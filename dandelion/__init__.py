import argparse

__version__ = '0.7.0'

def print_separator(text, width=70):
    border = "╔" + "═" * (width-2) + "╗"
    
    total_symbols_len = width - len(text) - 4  
    half_len = total_symbols_len // 2
    left_symbol = "║" + " " * (half_len - 1)
    right_symbol = " " * (total_symbols_len - half_len - 1) + "║"
    separator = left_symbol + '  ' + text + '  ' + right_symbol
    
    end = "╚" + "═" * (width-2) + "╝"
    print("\n\n" + border)
    print(separator)
    print(end + "\n\n")

def merge_args_with_defaults(module_parser, custom_args):
    """
    Merge custom arguments with module defaults.
    Args:
    - module_parser: the module parser function
    - custom_args: dictionary of custom arguments

    Returns:
    - argparse.Namespace: merged namespace of arguments
    """
    
    parser = module_parser()
    for action in parser._actions:
        if action.required:
            action.required = False

    defaults = vars(parser.parse_args([]))
    defaults.update(custom_args)

    for action in parser._actions:
        if not action.required and action.dest in custom_args:
            action.required = True

    return argparse.Namespace(**defaults)