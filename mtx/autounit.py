#!/usr/bin/env python3

"""
recursive automatic jemjive script

Python 2/3 compatible (well, apart from the print statement)

This script starts in the current directory and recursively
enters subdirectories to initiate a (JemJive) compile build,
execute a (JemJive) binary, a sphinx html document creation build,
or perform a unittest on output files. A combination of these tasks is
also possible.

Build:

Searches for file 'Makefile' in directory 'src' and calls 'make opt',
which is specific for JemJive, but it can be modified for use with
other build environments.

Execute:

Searches for JemJive input files ending with .pro in directories
containing (sub)project data.

 In case it finds a .pro file it searches for the directory ../src
 and looks for a program binary with the extension '-opt'. If found,
 it executes the binary in the directory where the .pro file is found.
 This is done via a 'with' context switch.

 Note that it will do this for every .pro file it finds, so take care
 to separate data files if you have more than one .pro file in a
 project directory!

Document-build:

Searches for file 'Makefile' in directory 'sphinx' and calls 'make html'

Unittest:

searches for either a directory "unittest" or a file extension .unittest

 In case it finds the directory "unittest" it will make a list of all
 files with the same name in the directory unittest and the directory
 immediately above. It will then compare the files with numdiff with
 a relative error of 1e-6

 In case it finds files with extension ".unittest" it will look for
 files without that extension in the same directory and compare them
 with numdiff with a relative error of 1e-6

Copyright 06-09-2016 Frank Everdij F.P.X.Everdij@tudelft.nl

TODO:
- support for compressed files

"""

__version__ = "$Revision$"
# $Source$

import sys
import argparse
import os
from filecmp import dircmp
from subprocess import call


# ANSI escape sequences for simple color output
class bcolors:
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    OFF = '\033[0m'


class WorkingDirectory(object):
    def __init__(self, new_cwd):
        self.new_cwd = new_cwd
    
    def __enter__(self):
        self.old_cwd = os.getcwd()
        os.chdir(self.new_cwd)
    
    def __exit__(self, type, value, traceback):
        os.chdir(self.old_cwd)


is_exe = lambda fpath : os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def find_binary(root):
    while root:
        for binroot, dirs, files in os.walk(os.path.join(root, 'src')):
            for filename in files:
                if (filename[-4:] == '-opt'):
                    binary = os.path.join(os.path.abspath(binroot), filename)
                    if is_exe(binary):
                        return binary
        root=os.path.dirname(root)


def file_float_cmp(root, fa, fb, err, boolabs):
    errortype = 'absolute' if boolabs else 'relative'
    relabs = '-a' if boolabs else '-r'
    print ('    unittest: In', root)
    with WorkingDirectory(root):
        FNULL = open(os.devnull, 'w')
        try:
            status = call(['numdiff', '-q', relabs, str(err), fa, fb],
                            stderr=FNULL)
        except OSError:
            exit('Error in call to numdiff. Check if numdiff is installed')
        finally:
            FNULL.close()
        if status != 0:
            print (bcolors.RED,'--- Files', fa, 'and', fb,
                   'differ by more than', err, errortype, 'error.',
                   bcolors.OFF)
        else:
            print (bcolors.GREEN,'+++ Files', fa, 'and', fb, 'match.',
                   bcolors.OFF)
        return True if status == 0 else False


def new_file_float_cmp(root, fa, fb, err, boolabs):
    if len(fa):
        errortype = 'absolute' if boolabs else 'relative'
        relabs = '-a' if boolabs else '-r'
        print ('    unittest: In', root)
        with WorkingDirectory(root):
            FNULL = open(os.devnull, 'w')
            try:
                for i in zip(fa,fb):
                    status = call(['numdiff', '-q', relabs, str(err), i[0],
                             i[1]], stderr=FNULL)
                    if status != 0:
                        print (bcolors.RED,'--- Files', i[0], 'and', i[1],
                               'differ by more than', err, errortype, 'error.',
                               bcolors.OFF)
                    else:
                        print (bcolors.GREEN,'+++ Files', i[0], 'and', i[1],
                               'match.', bcolors.OFF)
            except OSError:
                exit('Error in call to numdiff. Check if numdiff is installed')
            finally:
                FNULL.close()
    return


def recursive_make(directory, excludedir, base, makeoption, resultstring):
    for root, dirs, files in os.walk(directory):
        dirs[:] = [d for d in dirs if d not in excludedir]
        if os.path.basename(root) == base:
            with WorkingDirectory(root):
                if os.path.exists('Makefile'):
                    print_result(call(['make', makeoption]),
                                 resultstring+': Make '+makeoption)


def print_result(status, string):    
    s = 'failed' if status else 'succesful'
    m = bcolors.RED + '---' if status else bcolors.GREEN + '+++'
    print (m, string, s, bcolors.OFF)


def main():
    parser = argparse.ArgumentParser(
        description='Performs automatic recursive operations on a'
                    ' JemJive project tree',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--build', action = 'store_true',
                       help = 'perform program builds')
    parser.add_argument('-e', '--execute', action = 'store_true',
                       help = 'execute programs with -opt extension')
    parser.add_argument('-d', '--docbuild', action = 'store_true',
                       help = 'perform sphinx document builds')
    parser.add_argument('-u', '--unittest', action = 'store_true',
                       help = 'perform unittest on output files')
    parser.add_argument('-o', '--onlyinunit', action = 'store_true', help =
                        'only build/execute when unittest files are present')
    parser.add_argument('-a', '--abserror', action = 'store_true', help =
                        'use absolute error instead of relative in unittest')
    parser.add_argument('errorval', metavar = 'error', type = float,
                        nargs = '?', default = 1e-6,
                        help = 'value for relative or absolute error')
    args = parser.parse_args(sys.argv[1:])

    if not (args.build | args.execute | args.docbuild | args.unittest):
        print (sys.argv[0], ": missing option(s). Use -h for help")

    # Make a set to omit version control directories from any walk
    excludedir = set(['.hg','.git'])


    # Do we need to build and execute only in unittestable projects?
    if args.onlyinunit:

        # Yes, so we make a list of projects directories by checking
        # if there is a 'src' directory directly beneath it
        jjproject=set([])
        for root, dirs, files in os.walk(os.path.curdir):
            dirs[:] = [d for d in dirs if d not in excludedir]
            if 'src' in dirs:
                jjproject.add(root)

        # Then for each project directory, check which ones have unittests
        jjprojectunit=set([])
        for j in jjproject:
            for root, dirs, files in os.walk(j):
                if os.path.basename(root) == 'unittest':
                    jjprojectunit.add(j)
                for name in files:
                    if os.path.splitext(name)[1] == '.unittest':
                        jjprojectunit.add(j)
                        
        # Then for each data directory, check which ones have unittests
        jjdataunit=set([])
        for j in jjprojectunit:
            for root, dirs, files in os.walk(j):
                for name in files:
                    if os.path.splitext(name)[1] == '.pro':
                        if 'unittest' in dirs:
                            jjdataunit.add(root)
                        for name2 in files:
                            if os.path.splitext(name2)[1] == '.unittest':
                                jjdataunit.add(root)

    else:

        # No, use the current path for all searches then
        jjprojectunit=set(os.path.curdir)
        jjdataunit=set(os.path.curdir)


    # Since we can combine option switches for all tasks we start with
    # build, then execute, then docbuild and finally unittest

    if args.build:
        for j in jjprojectunit:
            recursive_make(j, excludedir, 'src', 'opt', 'build')

    if args.execute:
        for j in jjdataunit:
            for root, dirs, files in os.walk(j):
                dirs[:] = [d for d in dirs if d not in excludedir]
                for name in files:
                    if os.path.splitext(name)[1] == '.pro':
                        binary = find_binary(root)
                        if binary:
                            with WorkingDirectory(root):
                                print_result(call([binary, name]),
                                             ' '.join(['execute:',
                                             os.path.basename(binary), name]))

    if args.docbuild:
        recursive_make(os.path.curdir, excludedir, 'sphinx', 'html', 'docbuild')

    filea = []
    fileb = []
    if args.unittest:
        for root, dirs, files in os.walk(os.path.curdir):
            dirs[:] = [d for d in dirs if d not in excludedir]
            filea[:] = []
            fileb[:] = []
            if 'unittest' in dirs:
                dcmp = dircmp(root, os.path.join(root, 'unittest'))
                filea = dcmp.common_files
                fileb = [os.path.join('unittest',n) for n in filea]
            for name in files:
                filename, extension = os.path.splitext(name)
                if extension == '.unittest':
                    if os.path.exists(os.path.join(root, filename)):
                        filea.append(filename)
                        fileb.append(name)   
            new_file_float_cmp(root, filea, fileb, args.errorval,
                               args.abserror)


if __name__ == '__main__':
    main()

            
