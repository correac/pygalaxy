import argparse


class ArgumentParser(object):
    """
    Class for handeling arguments that are passed to the script
    """

    # Input redshift.
    redshift: float
    # Option of calculation to make.
    output_option: str
    # Directory to output the data and figures to
    output_directory: str

    def __init__(self):

        parser = argparse.ArgumentParser(
            description="""General argument parser for pygalaxy pipeline."""
        )

        parser.add_argument(
            "-z",
            "--redshift",
            help="Input redshift. Example: 0",
            type=float,
            required=False,
        )

        parser.add_argument(
            "--calculatemcrit",
            help=(
                "Calculate critical halo mass. yes/no. Default no."
            ),
            type=str,
            required=False,
        )

        parser.add_argument(
            "-o",
            "--output-directory",
            help="Output directory where data is stored.",
            type=str,
            required=True,
        )

        args = parser.parse_args()
        if args.calculatemcrit is None:
            args.calculatemcrit = "no"

        self.redshift = args.redshift
        self.output_directory = args.output_directory
        self.output_option = args.calculatemcrit

        print("Parsed arguments:")
        print("---------------------\n")
        print(f"Redshift: {self.redshift}")
        print(f"Output directory: {self.output_directory}")
        print(f"Output option: {self.output_option}")
        print("")

