#!/usr/bin/env python
"""Create a new protocol with initialised main.py and README.md"""
from pathlib import Path

def main(name : str):
    """Initialises new folder and creates main.py and README.md
    """
    def main_string() -> str:
        """Creates the string of a main.py file
        """
        string = []
        string.append('#!/usr/bin/env python')
        string.append(f'"""Update description for {name}"""')
        string.append(f'\ndef main(dummy : None = None):\n    print("Automatically Generated Program for {name}")\n    return')
        string.append('\nif __name__ == "__main__":')
        string.append('    from argparse import ArgumentParser')
        string.append('    parser = ArgumentParser()')
        string.append('    parser.add_argument("--dummy", type=None, required=False, default=None)')
        string.append('    main(**vars(parser.parse_args()))')
        return '\n'.join(string)

    def readme_string() -> str:
        """Creates the string of a README.md file
        """
        string = []
        string.append(f"# {name.replace('-', ' ').title()}")
        string.append(f"## Summary")
        string.append(f"Enter the summary of the program here.")
        string.append(f"## Usage")
        string.append(f"```python main.py```")
        string.append(f"## Description")
        string.append(f"Create an extended description of the program here.")
        return '\n\n'.join(string)

    if len(name.split(" ")) != 1:
        raise NameError(f"Name is currently {name} but cannot contain spaces!")

    # create folder
    Path(f'./{name}').mkdir()

    # create main.py
    with open(f'{name}/main.py', 'w') as f:
        f.write(main_string())

    # create README.md
    with open(f'{name}/README.md', 'w') as f:
        f.write(readme_string())

    print(f"Successfullly create {name}/README.md and {name}/main.py!")
    return

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('name', type=str)
    main(**vars(parser.parse_args()))