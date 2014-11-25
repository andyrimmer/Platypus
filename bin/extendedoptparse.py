"""
Small extension of the optparse module to add required arguments. This avoids quite a lot
of extra argument-checking code.

This code is based on an example from the following url:

http://code.activestate.com/recipes/573441-extended-optparse-to-allow-definition-of-required-/
"""

import optparse
import glob
from copy import copy

###################################################################################################

def checkList(option, opt, value):
    """
    Utility function to parse a list from a string, or return an error
    """
    try:

        # Remove trailing white space
        value = value.strip(" ")

        if " " in value:
            raise OptionValueError, "Option list cannot have spaces. You must supply a comma-separated list."

        if "*" in value:
            return glob.glob(value)
        elif "," in value:
            return value.replace(" ", "").split(",")
        else:
            return value.split()

    except ValueError:
        raise OptionValueError("option %s: invalid comma-delimitedstring value: %r" % (opt, value))

###################################################################################################

class ExtendedOption(optparse.Option):
    """
    Extend the optparse.Option class by adding a 'required'
    attribute.
    """
    ATTRS = optparse.Option.ATTRS + ["required"]
    TYPES = optparse.Option.TYPES + ("list",)
    TYPE_CHECKER = copy(optparse.Option.TYPE_CHECKER)
    TYPE_CHECKER["list"] = checkList

    def __init__(self, *opts, **attrs):
        """
        Initialise optparse.Option attributes, and pre-pend the
        help string with (Required) if appropriate.
        """
        if attrs.get("required"):
            attrs['help'] = '(Required) ' + attrs.get('help', "")

        optparse.Option.__init__(self, *opts, **attrs)

###################################################################################################

class OptionParser(optparse.OptionParser):
    """
    Extend the optparse.OptionParser class
    """
    def __init__(self, **kwargs):
        """
        Initialise the optparse.OptionParser attributes, and set the
        class used for storing options to the ExtendedOption
        class.
        """
        kwargs['option_class'] = ExtendedOption
        optparse.OptionParser.__init__(self, **kwargs)

    def check_values(self, values, args):
        """
        Ensure that all required options are present, otherwise raise an
        exception.
        """
        for option in self.option_list:
            if option.required and not getattr(values, option.dest):
                self.error("Option %s is required\n\n" % (str(option)))

        return optparse.OptionParser.check_values(self, values, args)

###################################################################################################
