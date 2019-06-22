"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

import json
import os

from pymongo.errors import BulkWriteError, DuplicateKeyError

from macrocycles.utils.database import MongoDataBase
from macrocycles.utils.exceptions import (SavingMongoError, SavingSQLError,
                                          WritingJsonError, WritingTxtError)
from macrocycles.utils.utils import DATA_DIR, PROJECT_DIR, Flags


class Base():
    """
    Class from which all other classes in package will inherit from. Handles output of result data.

    Attributes:
        project_dir: The filepath to the root project directory.
        data_dir: The filepath to the data directory.
        fp_out: The filepath to the output file.
        mongo_db: A connection to the MongoDB where the result_data will be stored.
        collection: The collection in the MongoDB that result_data will be inserted into.
        sql_db: A connection to the SQL database where the result_data will be stored.
        logger: The logger of the child classes' module.
        result_data: A list to contain the result_data.
        flags: A named tuple containing the following:
            json_flag: If true, save_data() writes result_data to a json file.
            txt_flag: If true, save_data() writes result_data to a txt file.
            mongo_flag: If true, save_data() writes result_data to the Mongo database.
            sql_flag: If true, save_data() writes result_data to the SQL database.
    """

    def __init__(self, f_out, collection, sql_db, logger, flags=Flags(False, False, False, False)):

        # I/O
        self.project_dir = PROJECT_DIR
        self.data_dir = DATA_DIR
        self.fp_out = os.path.join(DATA_DIR, f_out)
        self.mongo_db = MongoDataBase(logger=logger)
        self.collection = collection
        self.sql_db = sql_db
        self.logger = logger

        # data
        self.result_data = []

        # flags
        self.flags = flags

    def save_data(self):
        """
        Top level function for saving the data stored in self.result_data to all specified data formats. Calls helper
        function self.save_all_formats().

        Returns:
            bool: True if successful.
        """

        try:
            self.save_all_formats()
        except (WritingJsonError, WritingTxtError, SavingMongoError, SavingSQLError):
            pass
        else:
            if True not in self.flags:
                self.logger.warning('No data saved. Please set a data flag specifying which data format to save to.')
            return True

        return False

    def save_all_formats(self):
        """
        Helper function of self.save_data(). Sequentially calls the following functions: self.write_json(),
        self.write_txt(), self.save_mongo_db(), and self.save_sql_db(), in order to specifically diagnose any errors
        raised during exectution of each of these functions.

        Raises:
            WritingJsonError: Failure to write data in json format.
            WritingTxtError: Failure to write data in txt format.
            SavingMongoError: Failure to save data to the Mongo database.
            SavingSQLError: Failure to save data to the SQL database.
        """

        if self.flags.json_flag:
            try:
                self.write_json()
            except TypeError:
                self.logger.exception('Failed to write result_data to a json file')
                raise WritingJsonError
            else:
                self.logger.info(f'Successfully wrote data to {self.fp_out.split("/")[-1]}.json!')

        if self.flags.txt_flag:
            try:
                self.write_txt()
            except TypeError:
                self.logger.exception('Failed to write result_data to a txt file')
                raise WritingTxtError
            else:
                self.logger.info(f'Successfully wrote data to {self.fp_out.split("/")[-1]}.txt!')

        if self.flags.mongo_flag:
            try:
                self.save_mongo_db()
            except (DuplicateKeyError, ValueError):
                self.logger.exception('Failed to save result_data to the Mongo database')
                raise SavingMongoError
            except BulkWriteError as err:
                self.logger.exception(f'Failed to save result_data to the Mongo database\n{err.details}')
                raise SavingMongoError
            else:
                self.logger.info(f'Successfully saved data to the Mongo database!')

        if self.flags.sql_flag:
            try:
                self.save_sql_db()
            except Exception:
                self.logger.exception('Failed to save result_data to the SQL database')
                raise SavingSQLError
            else:
                self.logger.info(f'Successfully saved data to the SQL database!')

    def write_json(self):
        """
        Writes data stored in self.result_data to a json file specified by self.fp_out.
        """

        with open(self.fp_out + '.json', 'w') as f:
            json.dump(self.result_data, f)

    def write_txt(self):
        """
        Writes data stored in self.result_data to a txt file specified by self.fp_out.
        """

        with open(self.fp_out + '.txt', 'w') as f:

            # write data properties
            f.write(','.join(self.result_data[0].keys()) + '\n')

            # write data
            for data_point in self.result_data:
                data_string = [str(field) for field in data_point.values()]
                f.write(','.join(data_string) + '\n')

    def save_mongo_db(self):
        """
        Saves data stored in self.result_data to the Mongo database specified by self.mongo_db.

        Returns:
            bool: True if successful.
        """

        return self.mongo_db.insert(self.collection, self.result_data)

    def save_sql_db(self):
        """
        Saves data stored in self.result_data to the SQL database specified by self.sql_db.

        Returns:
            bool: True if successful.
        """

        return self.sql_db.insert(self.result_data)
