# compatibility.py
from inspect import signature

def formatargspec(*args, **kwargs):
    sig = signature(args[0])
    return str(sig)