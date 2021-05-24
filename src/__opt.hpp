#pragma once




#include <iterator>     // std::distance
#include <sstream>      // std::istringstream
#include <stdexcept>    // std::runtime_error
#include <string>




namespace nsl
{ // namespace nsl




/**
*******************************************************************************
** class __base_opt
*******************************************************************************
**/




class __base_opt
{ // class __base_opt


protected: /******************************************************************/


    std::string const& _prog_name;
    std::string _opt_name;

    char** _beg;


public: /*********************************************************************/


    /**
    *****************************************************************
    ** constructors
    *****************************************************************
    **/

    // Most complete constructor
    __base_opt( std::string const& prog_name,
                std::string const& opt_name ,
                char** beg
                );


    // If the option value is not specified
    __base_opt( std::string const& prog_name,
                std::string const& opt_name
                );


    /**
    *****************************************************************
    ** operators
    *****************************************************************
    **/

    operator bool() const;


    /**
    *****************************************************************
    ** member functions
    *****************************************************************
    **/

    auto is_specified() const
        -> bool;


    /**
        Function "convert" which simplifies the interface of the stream conversion
        with "std::istringstream" so that a failed conversion makes the function
        throws an exception. This only enables a conversion from anything bindable
        with a "std::string" to a built-in type: the first argument has thus to be
        bindable to a "std::string".

        The second argument of the function is for debugging: the exception message
        will take into account the name of the variable to which the conversion
        failed.

            @dev: using metaprogramming and SFINAE to partially specialise the
            template function.

    **/

    // implicit conversion (string to string)
    template < typename To, typename From >
    auto convert( From const& val ) const
        -> typename std::enable_if< std::is_same< To, std::string >::value,
                                    To /* real return type */ >::type
    {
        return val;
    }


    // boolean conversion (string to bool)
    template < typename To, typename From >
    auto convert( From val ) const
        -> typename std::enable_if< !std::is_same< To, std::string >::value &&
                                    std::is_same< To, bool >::value,
                                    To /* real return type */ >::type
    {
        // Converting to a string to ease string comparisons ////////////
        std::string tmp( std::move(val) );

        // Special checking /////////////////////////////////////////////
        if( tmp == "true" )
            return true;
        else if( tmp == "false" )
            return false;

        // Default conversion ///////////////////////////////////////////
        else
            return __convert< To >( tmp );
    }

    // most general conversion (string to any type)
    template < typename To, typename From >
    auto convert( From const& val ) const
        -> typename std::enable_if< !std::is_same< To, std::string >::value &&
                                    !std::is_same< To, bool >::value,
                                    To /* real return type */ >::type
    {
        return __convert< To >( val );
    }



protected: /******************************************************************/

    /**
    *****************************************************************
    ** private functions
    *****************************************************************
    **/

    auto __throw_if_unspecified() const
        -> void;

    auto __throw_if_invalid_conversion( std::istringstream const& str ) const
        -> void;

    auto __throw( std::string const& what ) const
        -> void;


    template < typename ReturnType >
    auto __convert( std::string const& value ) const
        -> ReturnType
    {
        // DONE: benchmark shows it against "no static and no clear" //////////
        static std::istringstream is;

        // Initialising /////////////////////////////////////////////////
        is.str( value );
        auto res = ReturnType{};

        // Converting ///////////////////////////////////////////////////
        is >> res;

        // Throwing if the conversion failed ////////////////////////////
        __throw_if_invalid_conversion( is );

        // TODO: benchmark it against "no static and no clear" //////////
        is.clear();

        // NRVO /////////////////////////////////////////////////////////
        return res;
    }


}; // class __base_opt *******************************************************/




/**
*******************************************************************************
** class __unnamed_opt
*******************************************************************************
**/



class __unnamed_opt: public __base_opt
{ // class __unnamed_opt


public: /*********************************************************************/


    /**
    *****************************************************************
    ** constructors
    *****************************************************************
    **/

    // Most complete constructor
    __unnamed_opt( std::string const& prog_name,
                   std::string const& opt_name ,
                   char** beg
                   );


    // If the option value is not specified
    __unnamed_opt( std::string const& prog_name,
                   std::string const& opt_name
                   );


    /**
    *****************************************************************
    ** operators
    *****************************************************************
    **/

    // Conversion operators -------------------------------
    template < typename T >
    operator T() const
    {
        // An unnamed option is mandatorily specified; otherwise an error would
        // have been thrown before.
        //// __throw_if_unspecified();

        return convert< T >( *_beg );
    }


}; // class __unnamed_opt ****************************************************/




/**
*******************************************************************************
** class __opt
*******************************************************************************
**/




// Forward declaration for %__check_bindable
class __opt;


/**
***********************************************************
** class __check_bindable (helper class)
***********************************************************
**/
template < typename T >
class __check_bindable
{ // class __check_bindable


private: /************************************************/

    T _val;
    __opt const& _opt;

public: /*************************************************/

    __check_bindable( __opt const& opt, T val ):
        _val( std::move(val) ), _opt( opt )
    {
        /* No further construction. */
    }


    template < typename U >
    operator U() const;

}; // class __check_bindable *****************************/




class __opt: public __base_opt
{ // class __opt


private: /********************************************************************/


    char** _end;


public: /*********************************************************************/


    /**
    *****************************************************************
    ** constructors
    *****************************************************************
    **/

    // Most complete constructor
    __opt( std::string const& prog_name,
           std::string const& opt_name ,
           char** beg                  ,
           char** end
           );


    // If the option value is not specified
    __opt( std::string const& prog_name,
           std::string const& opt_name
           );


    /**
    *****************************************************************
    ** operators
    *****************************************************************
    **/

    // Conversion operators -------------------------------
    template < typename T >
    operator T() const
    {
        __throw_if_unspecified();
        __throw_if_missing_values();
        __throw_if_too_many_values();

        return convert< T >( *_beg );
    }

    template < typename T,
               template <typename,typename> class Container
               >
    operator Container< T, std::allocator< T > >() const
    {
        __throw_if_unspecified();
        __throw_if_missing_values();

        return convert_n< T, Container >( _beg, _end );
    }

    // Default value operator -----------------------------
    template < typename T >
    auto operator| ( T rhs ) const
    {
        /// @dev : we want to use %operator T(), and not to convert to a rvalue
        ///        which could be of any type. Hence the helper class.
        return __check_bindable< T >( *this, rhs );
    }


    /**
    *****************************************************************
    ** member functions
    *****************************************************************
    **/

    /**
        The container version: it converts every range of element defined by the
        iterators passed in the arguments. Every element has to be convertible

            @dev: the order of the template arguments is important: the first must
            be as the other overloads.

        TODO: - make it as std functions? (i.e. not allocating) -> not possible
                to use auto, but it separates interfaces...

    **/
    template < typename T,
               template <typename,typename> class Container,
               typename It // deduced from function call
               >
    auto convert_n( It beg, It end ) const
        -> Container< T, std::allocator<T> >
    {
        // Initialising /////////////////////////////////////////////////
        auto res = Container< T, std::allocator<T> >{};
        res.reserve( std::distance( beg, end ) );

        // Converting each element //////////////////////////////////////
        for( auto it = beg; it != end; ++it )
            res.push_back( convert< T >( *it ) );

        // NRVO /////////////////////////////////////////////////////////
        return res;
    }


private: /********************************************************************/


    /**
    *****************************************************************
    ** private functions
    *****************************************************************
    **/

    auto __throw_if_missing_values() const
        -> void;

    auto __throw_if_too_many_values() const
        -> void;


}; // class __opt ************************************************************/




// %operator T() of %__check_bindable.
template < typename From >
template < typename To >
__check_bindable< From >::operator To() const
{
    // This is an error which is not depending on the options of the program but an error of the program.
    static_assert( std::is_convertible< From, To >::value, "invalid default value" );

    if( _opt.is_specified() )
        return _opt.operator To();
    else
        return _val;
}


// %operator| to reverse the order of default value and options.
template < typename T >
auto operator| ( T lhs, __opt const& rhs )
{
    return rhs | lhs;
}




} // namespace nsl


