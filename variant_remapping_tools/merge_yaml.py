#! /usr/bin/env python3
import argparse

import yaml


def merge_two_dict(d1, d2):
    result = {}
    for key in set(d1) | set(d2):
        if isinstance(d1.get(key), dict) or isinstance(d2.get(key), dict):
            result[key] = merge_two_dict(d1.get(key, dict()), d2.get(key, dict()))
        else:
            result[key] = d1.get(key, 0) + d2.get(key, 0)
    return result


def merge_yaml_files(input_yamls, output_yaml):
    output = {}
    for input_yaml in input_yamls:
        with open(input_yaml) as open_input:
            data = yaml.safe_load(open_input) or {}
            output = merge_two_dict(output, data)

    with open(output_yaml, 'w') as open_output:
        yaml.safe_dump(output, open_output)


def main():
    description = ('Merge multiple yaml file containing stats by summing the overlapping fields')

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--inputs', type=str, required=True, nargs='+',
                        help='YAML files containing the input summary metrics')
    parser.add_argument('--output', type=str, required=True,
                        help='YAML files containing the output summary metrics')
    args = parser.parse_args()

    merge_yaml_files(args.inputs, args.output)


if __name__ == '__main__':
    main()
