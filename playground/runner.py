def main():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Creates a unique set of molecules by attaching varying length alkyl '
                                     'chains to all elgible positions on the parent side chain. Alkyl chains include '
                                     'methyl, ethyl, and propyl.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--f_in', dest='f_in', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file relative to default input directory defined in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='side_chains',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')

    args = parser.parse_args()

    # check for proper file specifications
    if args.input in ['json', 'txt']:
        extension = args.f_in.split('.')[-1]
        if args.input != extension:
            LOGGER.error('File extension of the input file does not match the specified format')
            raise OSError('File extension of the input file does not match the specified format')

    # configure I/O
    input_flags, output_flags = set_flags(args.input, args.output)
    f_in = [args.f_in] if args.input in ['json', 'txt'] else ['']

    # create class and perform operations
    modifier = SideChainModifier(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                 output_flags=output_flags, no_db=args.no_db)
    if modifier.load_data() and modifier.diversify():
        return modifier.save_data()

    return False


if __name__ == '__main__':
    main()
