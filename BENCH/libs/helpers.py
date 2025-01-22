##########################################################
#  CROCO psyclone build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
"""Provide some helper function for example to run commands and patch files."""

##########################################################
import re
import os
import sys
import shutil
import subprocess
import tempfile
import traceback
from termcolor import cprint
from timeit import timeit
from .messaging import Messaging
from contextlib import contextmanager


##########################################################
def display_run_error(command: str, outpout: str) -> None:
    dir = os.getcwd()
    keep_n_last_lines = 32
    truncated_output = "\n".join(outpout.split("\n")[-keep_n_last_lines:])
    cprint("-----------------------------------------------", "red")
    cprint(truncated_output, "dark_grey")
    cprint("-----------------------------------------------", "red")
    cprint(f"Error from command : {command}", "red")
    cprint(f"Error from workdir : {dir}", "red")
    cprint("-----------------------------------------------", "red")
    print("")


##########################################################
def print_exception(exception: Exception) -> None:
    dir = os.getcwd()
    cprint("-----------------------------------------------", "red")
    cprint("".join(traceback.format_exception(exception)), "dark_grey")
    cprint("-----------------------------------------------", "red")
    cprint(f"Error from command : {' '.join(sys.argv)}", "red")
    cprint(f"Error from workdir : {dir}", "red")
    cprint("-----------------------------------------------", "red")
    cprint(str(exception), "red")
    cprint("-----------------------------------------------", "red")
    print("")


##########################################################
def run_shell_command(command, logfilename=None, capture=True, show_on_error=True):
    """Print the command and run it. On failure it prints the output."""

    # print command
    Messaging.command(command)

    # capture or not capture
    if capture:
        # run
        result = subprocess.run(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if logfilename:
            with open(logfilename, "w") as log_file:
                log_file.write(result.stdout)
        # on failure
        if result.returncode != 0:
            if show_on_error:
                display_run_error(command, result.stdout)
            raise Exception("Fail to run command !")
    else:
        print("-----------------------------------------------")
        subprocess.run(command, shell=True, check=True)
        print("-----------------------------------------------")


def run_shell_command_time(command, logfilename=None, verbose: bool = False):
    """Print the command and run it by measuring time. On failure it prints the output."""
    # capture or not capture

    if verbose:
        return timeit(
            stmt=f"subprocess.run('{command}', check=True, shell=True)",
            setup="import subprocess",
            number=1,
        )
    else:
        if logfilename:
            log_fp = logfilename
        else:
            log_fp = tempfile.NamedTemporaryFile("w+").name

        try:
            return timeit(
                stmt=f"subprocess.run('({command}) 2>&1 > {log_fp}', check=True, shell=True)",
                setup="import subprocess",
                number=1,
            )
        except:
            display_run_error(command, log_fp)
            raise Exception("Fail to run command !")


def replace_in_file(in_path: str, out_path: str, pattern: str, replace_by: str):
    """Patch a file by replacing some patterns inside."""

    # print progress
    fname = os.path.basename(out_path)
    Messaging.step(f"Patching {fname}")

    # open & load
    with open(in_path, "r") as fp:
        content = fp.read()

    # patch
    content = content.replace(pattern, replace_by)

    # rewrite
    with open(out_path, "w+") as fp:
        fp.write(content)


def contains_one_of(seach_in: str, elements) -> bool:
    """check is the string contains one of the given elements."""

    # if search a simple string
    if isinstance(elements, str):
        return elements in seach_in
    elif isinstance(elements, list):
        for element in elements:
            if element in seach_in:
                return True
        return False
    else:
        raise Exception("Unsupported type !")


def normalize_line(line):
    """Normalize a line by removing extra spaces and trimming."""
    if line is None:
        return ""
    return re.sub(r"\s+", " ", line).strip()


def apply_rule(fname, lines, rule):
    """
    Apply a single rule to a list of lines in memory.

    Parameters:
    - lines: a list of strings (lines from the file)
    - rule: a dictionary with keys:
        - "what" (str, required): The line to match (after normalization) for certain operations.
        - "by" (str): Line used in 'replace' mode.
        - "insert" (str or list): Line(s) to insert in 'insert-before', 'insert-after', 'insert-after-before' modes.
        - "mode" (str, required): One of 'replace', 'insert-before', 'insert-after', 'insert-after-before', 'insert-at-begin', 'insert-at-end'.
        - "descr" (str, optional): Description of the rule.
        - "before" (str, required for 'insert-after-before'): The line before which to insert.
        - "after" (str, required for 'insert-after-before'): The line after which to insert.
    """

    mode = rule.get("mode")
    what = normalize_line(rule.get("what", ""))
    descr = rule.get("descr", "No description provided")
    Messaging.step(f"Patching {fname} [ {descr} ]")

    # Determine the new line(s) to use based on the mode
    if mode == "replace":
        new_line = rule.get("by")
    else:
        new_line = rule.get("insert")

    # Basic validations
    if not mode:
        print(f"Error: No mode specified for rule '{descr}'. Skipping rule.")
        return lines
    if mode not in ["replace", "insert-before", "insert-after", "insert-after-before", "insert-at-begin", "insert-at-end"]:
        print(f"Error: Unknown mode '{mode}' in rule '{descr}'. Skipping rule.")
        return lines
    if not what and mode not in ["insert-at-begin", "insert-at-end"]:
        print(f"Error: 'what' is required for mode '{mode}' in rule '{descr}'. Skipping rule.")
        return lines
    if not new_line and mode != "replace":
        print(f"Error: 'insert' is required for mode '{mode}' in rule '{descr}'. Skipping rule.")
        return lines

    # Specific check for insert-after-before mode
    if mode == "insert-after-before":
        line_before = normalize_line(rule.get("before", ""))
        line_after = normalize_line(rule.get("after", ""))
        if not line_before or not line_after or not new_line:
            print(f"Error: 'before', 'after', and 'insert' must be provided for 'insert-after-before' in rule '{descr}'. Skipping rule.")
            return lines

    modified_lines = []
    inserted = False

    if mode == "replace":
        # Replace the line matching 'what' with 'by'
        for line in lines:
            if normalize_line(line) == what:
                if isinstance(new_line, list):
                    modified_lines.extend([nl + "\n" for nl in new_line])
                else:
                    modified_lines.append(new_line + "\n")
                inserted = True
            else:
                modified_lines.append(line)
        if not inserted:
            print(f"Warning: No line matched '{what}' for replacement in rule '{descr}'.")

    elif mode == "insert-before":
        # Insert a new line(s) before the line that matches 'what'
        for line in lines:
            if normalize_line(line) == what and not inserted:
                if isinstance(new_line, list):
                    modified_lines.extend([nl + "\n" for nl in new_line])
                else:
                    modified_lines.append(new_line + "\n")
                modified_lines.append(line)
                inserted = True
            else:
                modified_lines.append(line)
        if not inserted:
            print(f"Warning: No line matched '{what}' to insert before in rule '{descr}'.")

    elif mode == "insert-after":
        # Insert a new line(s) after the line that matches 'what'
        for line in lines:
            if normalize_line(line) == what and not inserted:
                modified_lines.append(line)
                if isinstance(new_line, list):
                    modified_lines.extend([nl + "\n" for nl in new_line])
                else:
                    modified_lines.append(new_line + "\n")
                inserted = True
            else:
                modified_lines.append(line)
        if not inserted:
            print(f"Warning: No line matched '{what}' to insert after in rule '{descr}'.")

    elif mode == "insert-after-before":
        # Insert a line(s) between two known lines: after 'after' and before 'before'
        line_before = normalize_line(rule.get("before", ""))
        line_after = normalize_line(rule.get("after", ""))

        found_after = False
        for i, line in enumerate(lines):
            current_norm = normalize_line(line)
            modified_lines.append(line)

            if found_after and current_norm == line_before and not inserted:
                if isinstance(new_line, list):
                    modified_lines[-1:-1] = [nl + "\n" for nl in new_line]
                else:
                    modified_lines.insert(len(modified_lines) - 1, new_line + "\n")
                inserted = True

            if current_norm == line_after:
                found_after = True

        if not inserted:
            print(f"Warning: Could not find the sequence 'after' then 'before' for rule '{descr}'. No insertion performed.")

    elif mode == "insert-at-begin":
        # Insert a new line(s) at the beginning of the file
        if not inserted:
            if isinstance(new_line, list):
                modified_lines.extend([nl + "\n" for nl in new_line])
            else:
                modified_lines.append(new_line + "\n")
            inserted = True
        modified_lines.extend(lines)

    elif mode == "insert-at-end":
        # Insert a new line(s) at the end of the file
        modified_lines.extend(lines)
        if not inserted:
            if isinstance(new_line, list):
                modified_lines.extend([nl + "\n" for nl in new_line])
            else:
                modified_lines.append(new_line + "\n")
            inserted = True

    return modified_lines



def patch_lines(file_path, rules):
    """
    Modify a file according to a list of rules.

    Parameters:
    - file_path: Path to the file to be modified.
    - rules: A list of dictionaries, each representing a rule.
    """
    # Backup the original file
    backup_path = f"{file_path}.backup"
    try:
        shutil.copyfile(file_path, backup_path)
        print(f"Backup created: {backup_path}")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found. Aborting.")
        return
    except Exception as e:
        print(f"Error while creating backup: {e}")
        return

    # Read the original file
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading the file: {e}")
        return

    # Apply each rule in-memory
    for rule in rules:
        lines = apply_rule(file_path, lines, rule)

    # Write the final result
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.writelines(lines)
    #        print(f"Modifications completed successfully for {file_path}.")
    except Exception as e:
        print(f"Error writing the file: {e}")


@contextmanager
def move_in_dir(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def apply_vars_in_str(value: str, vars: {}) -> str:
    # TODO: do not yet allow recursion but can ben done easly the way it is.

    # extract {case.value} patterns
    patterns = re.findall(r"\{[a-zA-Z0-9._-]+\}", value)

    # replace all
    if patterns:
        for pattern in patterns:
            # progress
            Messaging.step(f"Apply variable {pattern}...")

            # fetch value
            path = pattern.replace("{", "").replace("}", "").split(".")

            # start
            cursor = vars
            for element in path:
                cursor = cursor[element]

            # replace
            value = value.replace(pattern, cursor)

    # ok
    return value


def copy_tree_with_absolute_symlinks(src, dest):
    """
    Copies a directory tree, recreating folders and replacing files with
    absolute symbolic links.
    Existing symbolic links in the source are recreated as symbolic links in
    the destination.


    :param src: Path to the source directory tree.
    :param dest: Path to the destination directory tree.
    """
    for root, dirs, files in os.walk(src, followlinks=True):
        # Compute the corresponding path in the destination
        relative_path = os.path.relpath(root, src)
        dest_dir = os.path.join(dest, relative_path)

        # Recreate folders in the destination
        os.makedirs(dest_dir, exist_ok=True)

        # Process directories and files in the current folder
        for name in dirs + files:
            src_path = os.path.join(root, name)
            dest_path = os.path.join(dest_dir, name)

            # Check if the source is a symbolic link
            if os.path.islink(src_path):
                # Get the target of the symbolic link
                link_target = os.readlink(src_path)
                # Convert to absolute path if necessary
                link_target = os.path.join(os.path.dirname(src_path), link_target)
                link_target = os.path.abspath(link_target)

                # Check if the target is a file or directory
                if os.path.isfile(link_target):
                    # Create the symbolic link in the destination
                    try:
                        os.symlink(link_target, dest_path)
                    except FileExistsError:
                        print(
                            f"The symbolic link {dest_path} already exists, link not created."
                        )
                    except OSError as e:
                        print(f"Error creating symbolic link for {src_path}: {e}")

                elif os.path.isdir(link_target):
                    pass
            elif os.path.isfile(src_path):
                # Create a symbolic link for a regular file
                try:
                    os.symlink(os.path.abspath(src_path), dest_path)
                except FileExistsError:
                    print(f"The file {dest_path} already exists, link not created.")
                except OSError as e:
                    print(f"Error creating link for {src_path}: {e}")
            elif os.path.isdir(src_path):
                # Skip folders as they are already handled by os.walk
                pass
