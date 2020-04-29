import json


def text_load(filepath):
    with open(filepath, 'r') as file:
        return file.readlines()


def text_save(data, filepath):
    with open(filepath, 'w') as file:
        for line in data:
            file.write(line)


def json_load(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)


def json_save(data, filepath):
    with open(filepath, 'w') as file:
        json.dump(data, file)
