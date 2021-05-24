#pragma once




#include <stdexcept>




namespace nsl
{ // namespace basic_function




struct exit_exception : public std::runtime_error
{
    exit_exception() : std::runtime_error( "" ) {}
};




} // namespace basic_function
