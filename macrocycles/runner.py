import argparse


def side_chain_generator():
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


def monomer_generator():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Generates monomers by combining side chains and backbone together. If'
                                     ' the set of generated monomers are to be required in peptide generation, need to '
                                     'spcify option --required.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--bb_file', dest='bb_file', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing backbone data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--sc_file', dest='sc_file', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing side chain data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='custom',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')
    parser.add_argument('--required', dest='required', action='store_true',
                        help='Determines whether the monomer will be in the required pool for generating peptides. '
                        'Defaults to False.')

    args = parser.parse_args()

    # check for proper file specifications
    if args.input in ['json', 'txt']:
        extensions = [args.bb_file.split('.')[-1]]
        extensions.append(args.sc_file.split('.')[-1])
        if not all([args.input == extension for extension in extensions]):
            LOGGER.error('File extension of the input files does not match the specified format')
            raise OSError('File extension of the input files does not match the specified format')

    # configure I/O
    input_flags, output_flags = set_flags(args.input, args.output)
    f_in = [args.sc_file, args.bb_file] if args.input in ['json', 'txt'] else ['']

    generator = MonomerGenerator(required=args.required, f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                 output_flags=output_flags, no_db=args.no_db)
    if generator.load_data() and generator.generate_monomers():
        return generator.save_data()

    return False


def generate_peptides():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(
        description='Connects specified number of monomers to form a peptide. Takes input file(s) containing '
        'monomer SMILES strings and outputs to a json file the peptides as SMILES strings. Input and output are json '
        'files because this script will run on a super computer')
    parser.add_argument('length', type=int, help='The length of the peptide in monomers.')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--f_in', dest='f_in', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing backbone data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='custom',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')
    parser.add_argument('--num_peptides', type=int, dest='num_peptides', help='The number of peptides to generate')
    parser.add_argument('--num_jobs', dest='num_jobs', default=1, type=int,
                        help='The number of jobs to run on a job array.')
    parser.add_argument('--job_id', dest='job_id', default=1, type=int, help='The job number or ID.')

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

    generator = PeptideGenerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                 output_flags=output_flags, no_db=args.no_db)
    # t = time()
    if generator.load_data() and generator.generate_peptides(args.length, args.num_peptides,
                                                             args.num_jobs, args.job_id):
        # print(len(generator.result_data))
        # print(generator.result_data[0])
        # print(time() - t)
        return generator.save_data()

    return False


def tp_hybrid_generator():
    """
    Driver function. Parses arguments, constructs class, and performs operations on data.
    """

    parser = argparse.ArgumentParser(description='Connects each peptide defined in the input file to a template and '
                                     'write the resulting molecule to file as a SMILES string')
    parser.add_argument('input', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies the format that the input data is in.')
    parser.add_argument('output', nargs='+', choices=['json', 'txt', 'mongo', 'sql'],
                        help='Specifies what format to output the result data in.')
    parser.add_argument('--bb_file', dest='bb_file', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing backbone data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--sc_file', dest='sc_file', required='json' in sys.argv or 'txt' in sys.argv,
                        help='The input file containing side chain data relative to default input directory defined '
                        'in config.py.')
    parser.add_argument('--f_out', dest='f_out', default='custom',
                        help='The output file relative to the default output directory defined in config.py.')
    parser.add_argument('--no_db', dest='no_db', action='store_true',
                        help='Turns off default connection that is made to the database.')

    args = parser.parse_args()

    # check for proper file specifications
    if args.input in ['json', 'txt']:
        extensions = [args.bb_file.split('.')[-1]]
        extensions.append(args.sc_file.split('.')[-1])
        if not all([args.input == extension for extension in extensions]):
            LOGGER.error('File extension of the input files does not match the specified format')
            raise OSError('File extension of the input files does not match the specified format')

    # configure I/O
    input_flags, output_flags = set_flags(args.input, args.output)
    f_in = [args.sc_file, args.bb_file] if args.input in ['json', 'txt'] else ['']

    generator = TPHybridGenerator(f_in=f_in, f_out=args.f_out, input_flags=input_flags,
                                  output_flags=output_flags, no_db=args.no_db)
    if generator.load_data() and generator.generate_tp_hybrids():
        # print(len(generator.result_data))
        # print(generator.result_data[:3])
        return generator.save_data()

    return False


if __name__ == '__main__':
    main()
