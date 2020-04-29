import json
import csv


def load_text(filepath):
    with open(filepath, 'r') as file:
        return file.readlines()


def save_text(data, filepath):
    with open(filepath, 'w') as file:
        for line in data:
            file.write(line)


def load_json(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)


def save_json(data, filepath):
    with open(filepath, 'w') as file:
        json.dump(data, file)


def load_csv(filepath):
    with open(filepath, 'r') as file:
        return list(csv.reader(file, delimiter=','))
