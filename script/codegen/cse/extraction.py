from argparse import ArgumentParser
from source.cse_modules import extractDynamicsExpressions, extractOperatorExpressions

# Parse name of dynamics to be optimized
parser = ArgumentParser()
parser.add_argument('target', help='full C++ typename of the target instance')
parser.add_argument('template', help='path and filename to mako template')
parser.add_argument('input', help='path and filename to .txt file containing target info')
parser.add_argument('cpp_output', help='path where the .cpp files are created')
parser.add_argument('out_output', help='path where the .out files are created')
parser.add_argument('type', help='dynamics or operator')
args = parser.parse_args()

# Create .cpp files to extract expression tree
if args.type == 'dynamics':
    extractDynamicsExpressions(args.target, args.template, args.input, args.cpp_output, args.out_output)
elif args.type == 'operator':
    extractOperatorExpressions(args.target, args.template, args.input, args.cpp_output, args.out_output)
else:
    print("Invalid type.")
