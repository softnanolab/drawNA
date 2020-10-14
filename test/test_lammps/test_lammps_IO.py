#!/usr/bin/env python
import os

from drawNA.oxdna import System

ROOT = "/".join(os.path.abspath(__file__).split("/")[:-1])

def test_read_data():
    return

def test_read_dump():
    f_template = f'{ROOT}/data/dump.test.{{}}.lammpstrj'
    return

def test_write_data():
    return

if __name__=='__main__':
    test_read_data()
    test_read_dump()
    test_write_data()

