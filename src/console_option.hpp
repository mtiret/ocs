/******************************************************************************

NSL (Non Standard Library)

Console option treatment.

Last update : 22/04/2016
Beginning   : 30/06/2013

Author: Mathieu Tiret




*************************************************
** DESCRIPTION
*************************************************


It often occurs for console applications to read options before running in a
way the user wanted it to. This file helps reducing the redundant code doing
such a thing: it provides a core object, "console_option", with which one
function allows to extract the value(s) of the options.

The way to interpret the console options shapes the program and the way it is
used, therefore this module also aims to shape a command line standard. This
starts with observations: among the options, there are option names and option
values. There are only two possible option names for an option: the @short and
the @long version starting respectively with one or two hyphens.

There are also a special requirement about the @unnamed options (which are
mandatory): the user must put them before any of the @named options, namely as
the first options of the command line. Some already developed programs can read
these arguments from anywhere, but in a sense it makes the command line
difficult to read (since it is possible for some options to take a vector as
value). We choose here not to allow this kind of program to be developed, so
that the unnamed options are just after the program name in a command line.
Here is the standard this module tries to make the programs follow:

    @standard :
        `program` `unnamed options`... [option name][option value(s)]...
            `` : mandatory
            [] : optional

Finally, this module aims to shape the source code by reducing not only
redundancy between programs but also within the program which uses this module.
Indeed, all the functions and objects have been designed to make it possible
for the developers to write only once all the information concerning an option,
which therefore reduces code dependency between functions in the program. It is
thus recommended to call the functions with rvalues (it is possible to call
them with lvalues though).

-------------------------------------------------

 @summary :

    1. Reduce redundancy between programs for reading program options.

    2. Shape the program so that the command line follows a kind of standard
       program call.

    3. Reduce dependency between parts of the program, by treating every
       information about options and providing anything a developer could need.


*************************************************
** VOCABULARY
*************************************************


We will hereafter use the following terms:

- @dev: a programmer which uses this module.

- @user: a programmer which uses a program (in a command line) implemented with
  this module.

- @input: every bunch of characters which are separated by a space and which
  describes a command line (i.e. a program call). There are two types of
  inputs: the input names and the input values.

    - input name: a string which starts with either two hyphens (long version),
      or with one hyphen (short version).

    - input value(s): every string separated by a space character and between
      two input names (i.e. every input but the input name).

  There are several pattern for the inputs:

    - just an input name (hence for binary inputs)  : --input
    - input with value (white space between)        : --input 42
    - input with value (equal sign between)         : --input=42
    - input with values (white space between)       : --input 42 42 42
  Since it is more difficult to read, if an equal sign separates an input
  name and several values, it is considered as an error. @todo

- @options: they are sets of an input name and its value(s). An option could be
  mandatory or not (the unnamed options are mandatory). The developer specifies
  in its program if the option is mandatory or not, or its mandatoryness can
  follow a logical expression (more flexible).

- read an option: find the option name among the inputs and extract its value
  into a "std::string".

- specify an option: the user has effectively specified in the command line the
  option name and its value(s).


*************************************************
** LIBRARY USAGE
*************************************************


In order to make the interface simpler, there are only two functions the
dev should use to extract the inputs: "exit_if" and "read". Therefore, there
are several overloads of these functions.

Roughly, this file provides a smart mapping and an interface which goes with
it, so that developers only need the input name (as a key) to read its
value(s). It is noteworthy that the hyphen(s) are not removed from the keys.

Since the user should be using rvalues, most of the function arguments are
taken by value (copy elision makes it fast enough).
    @todo : use universal reference wherever it is possible.


-------------------------------------------------


- "exit_if":

        @prototype 1. auto exit_if( string s, string l, string doc ) -> void
        @prototype 2. auto exit_if( string s, string l, function   ) -> void

        where "s" and "l" are respectively the short and the long version of an
        option name (the order does not matter).

  This is a function which exits the program (by throwing an exception) if an
  option ("s" or "l") is specified. The behaviour before exiting the program
  depends on which overload is used. The former overload automatically
  generates a documentation, prints it in an error stream ("std::clog") and
  throws an exception; the latter calls a dev defined function (without any
  argument) and throws an exception.

  The exception thrown is of a special type ("nsl::exit_exception"), but it
  inherits from "std::exception": it is thus up to the eveloper to treat it
  differently than the other exceptions or not.

  Usually this function is used to exit the program if a help or a version
  option is specified.


-------------------------------------------------


- "read":

        @prototype 1. auto read( string uname, size_t n  ) -> T
        @prototype 2. auto read( string s, string l      ) -> T

    where T is the type of the return type (the dev does not have to deal
    with).

  This is a function which reads the options, stores the option values and
  converts them into a specified type. This function should be used as a rvalue
  for initialising a parameter. The option value can be missing if the
  parameter's type is a boolean: the value of a binary option is thus false if
  not specified and true if it is.
    @dev: this function do not convert anything, it returns a hidden class type
          "__opt" which converts the value depending on the return type
          ("operator T()") and check several other things.

  Each different overload reads different kind of options:

        1. reads the nth unnamed arguments, which internal name is "uname".
           (no conflict because there is not hyphen).
        2. reads a named option.

  Both of the overloads have versions which takes as third and/or fourth
  argument a documentation description and/or a tranformer (doc + transformer,
  doc only or transformer only). The description will be used to generate
  the documentation and the transformer makes it possible to customise the
  converting process by first applying the tranformer function to the option
  values. This transformer must take as an argument one "std::string" and
  return one "std::string".


-------------------------------------------------


Mandatoriness of the options:

  It is possible for the dev to compell a certain degree of requirement of the
  option.

    - the option is by default mandatory, so that if the option name is not in
      the command line an exception is thrown, except if the option is a binary
      option, in which case no exception is thrown.

    - it is possible to provide a default value by assigning at the end of the
      line the default value:

        int N = a.read( "-short", "--long" ) = 5;
        @todo int N = 5 | a.read( "-short", "--long" )

    - it is possible to have a logical relationship between the options, so
      that the mandatoriness is more flexible:

        int N = a.read( "s1", "l1" ) & a.read( "s2", "l2" );

      Each call to the function "read" can be interpreted as a boolean which is
      true if the option name is effectively found in the command line. If the
      boolean combination results to false, an exception is thrown.

      Here are the available boolean operators:
        & (and), | (exclusive or, xor), ^ (xor) and usage of the parentheses.

      The logical operator "or" does not have a sense since at the end the
      object has to return only one value. Hence, only "xor" is considered.
      In a case of ambiguity, the operator "^" is also overloaded as a "xor".

  The API compells the dev to specify the option type at the beginning of the
  line (and not using auto, except if a validator or a default value is
  provided).


-------------------------------------------------


  Before using all those functions, this object has to be initialised using the
  constructor which prototype is as follows:

        @prototype console_option( int argc, char** argv )


There are also two helper functions which could be useful in some case: "check"
and "get_documentation":

- "check" is a function which checks whether all the options in the command
  line have been read. It is useful if one wants to make its program stricter
  so that no computation is launched if the command line is not correct.

  If this function is called before some reading operations are finished, it
  will likely throw an exception for the not yet read options.

        @prototype void check()

- "get_documentation" returns the documentation description automatically
  generated by the function "read".

        @prototype auto get_documentation() -> string


-------------------------------------------------

 @summary of the possible usage:

    int N = a.read( "s1", "l1" ) | a.read( "s2", "l2" );

    auto N = static_cast< int >( a.read( "s1", "l1" ) | a.read( "s2", "l2" ) );

    auto N = a.read( "s1", "l1" ) | a.read( "s2", "l2" ) = 5;

    auto N = a.read( "l", "s", converter );


*************************************************
** EXCEPTIONS
*************************************************


There are several cases exceptions are thrown:

- if the value of a requested option is missing, except if the returned type
  is a boolean, namely a binary option.

- if an option is read several times.

- if an option is not read at all ( done by "check()" ).


*************************************************
** DEV
*************************************************


Actually, it is a mapping of ranges over "argv" with the option name as the
key.

"read" stores the unnamed options with a digit only as the key (which is the
position of the unnamed option in the command line), precedede by 'm' so that
there is not ambiguity with negative values.

If a rvalue version for a function is used, it works only if every parameter
are actually rvalues.


 @done :
    - 30/09/2014 : make it an object.
    - 11/02/2015 : "extract" should convert multiple value arguments into
                   vector of any type.
    - 28/02/2015 : core changes in the conception of the object
        - move all the converting work into an other object
        - no variadic functions
    - 04/03/2015 : adding several error checking in "_option_index()":
        - providing only "-" or "--" as an option name is an error.
        - providing a string starting with more than two hyphens as an option
          name is an error.
    - 04/03/2015 : why putting all the value in a string if it is afterwards
                   reconverted to a vector? Therefore change the data type to a
                   pair of vector of string and a boolean.
    - 25/03/2015 : working.
    - 24/07/2015 : redocumenting this file, and changing some standard.
    - 22/04/2016 : making the class easier: no combination...


 @todo :

    - consider several specification of an option as an error. Check if the
      different option names are specified and if it is the case, throw.

    - an automatic documentation function from the function extract.

    - make it possible to assign boolean with string ("true" and "false").

    - make a static array of message errors and add the program name, help...

    - make it possible to concatenate some short options

    - use regex or mapping?

    - use first lvalues const&, and if really needed use rvalue &&-> yes to
      reduce the number of redundancy

    - mutable and const?

    - change every iteration from 0 to an appropriate type

    - take care of the fact that argv is used without any copy of its values.
      Hence no move will be needed.

    - put more understandable error if __opt_conv fails (static_assert)

    - bool default value?

    - use the operator | for the default value (which should be specified first)



*******************************************************************************/




#pragma once




#include <unordered_map>
#include <vector>




#include "exception.hpp"
// includes: - std: stdexcept

#include "__opt.hpp"
// includes: - std: iterator, sstream, stdexcept, string




namespace nsl
{ // namespace nsl




class console_option
{ // class console_option


private: /********************************************************************/


    // Named option data type
    struct data_type
    {
        mutable bool used;  // has it been read/used?
        char **beg, **end;  // range of the values
    };

    // Unnamed option data type
    struct udata_type
    {
        mutable bool used;  // has it been read/used?
        char **val;         // value
    };

public: /*********************************************************************/


    /**
    *****************************************************************
    ** typedef
    *****************************************************************
        The key is the option name, the value is a tuple (a triplet) of the
        corresponding option value and a boolean indicating whether the option
        has been loaded.
    */
    using option_type = std::unordered_map< std::string, data_type >;



private: /********************************************************************/

    option_type _opt;
    std::vector< udata_type > _u_opt;

    std::string _doc;
    std::string _prog_name;


public: /*********************************************************************/


    /**
    *****************************************************************
    ** constructor
    *****************************************************************
        Reads the inputs of the program in the command line.
    **/


    console_option( int argc, char** argv );


    /**
    *****************************************************************
    ** exit_if
    *****************************************************************
        Performs the user defined function and throws an exception to end the
        program.

        @dev : using copy elision for "s" and "l".
    **/
    auto exit_if( std::string s,
                  std::string l,
                  std::string doc,
                  std::ostream& out
                  )
        -> void;


    template < typename F >
    auto exit_if( std::string s, std::string l, F fn )
        -> void
    {
        if( __is( s, l ) )
        {
            fn();
            throw nsl::exit_exception();
        }
    }


    /**
    *****************************************************************
    ** read
    *****************************************************************
        Reads a stored string from the map, and send it to a helper class
        "__opt".
    */

    /////////////////////////////////////////////////////////////////
    // Unnamed option overloads /////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    auto read( std::string uname    ,
               std::size_t n        ,
               std::string doc = ""
               )
        -> __unnamed_opt;

    /////////////////////////////////////////////////////////////////
    // Named option overloads ///////////////////////////////////////
    /////////////////////////////////////////////////////////////////


    auto read( std::string s    ,
               std::string l    ,
               std::string doc = ""
               )
        -> __opt;


    /**
    *****************************************************************
    ** check
    *****************************************************************
        Checks whether all the arguments read in the console have been loaded.
    */


    auto check() const
        -> void;


    /**
    *****************************************************************
    ** get_documentation
    *****************************************************************
        Returns the automatically generated documentation.
    */


    auto get_documentation() const
        -> std::string;


    auto command_line() const
        -> std::string;


private: /********************************************************************/


    /**
    *****************************************************************
    ** __is
    *****************************************************************
        Only eases the reading of the code: returns whether the option is in
        the command line.
    */


    auto __is( std::string const& s,
               std::string const& l
               ) const
        -> bool;


    /**
    *****************************************************************
    ** __find
    *****************************************************************
        Only eases the reading of the code (it is just a combination of std
        functions).
    */


    auto __find( std::string const& s,
                 std::string const& l
                 ) const
        -> auto
    {
        auto it = _opt.find( s );

        if( it == std::end( _opt ) )
            it = _opt.find( l );

        return it;
    }


    /**
    *****************************************************************
    ** __is_option_name
    *****************************************************************
        Returns whether the parameter is an option name. This functions checks also
        some inputs are ill formated.
    */


    inline
    auto __is_opt_name( char* opt ) const
        -> bool;


    /**
    *****************************************************************
    ** __throw
    *****************************************************************
        Only eases the reading of the code: it adds the program name in front
        of the exception message.
    */


    auto __throw_if_used( bool used, std::string const& name ) const
        -> void;


    auto __throw( std::string const& what ) const
        -> void;


    /**************************************************************************
    ***************************************************************************
    **************************************************************************/


}; // class console_option




} // namespace nsl












