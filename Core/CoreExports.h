#pragma once

#ifdef AEONTERRA_CORE_EXPORTS
    #define AEONTERRA_API __declspec(dllexport)
#else 
    #define AEONTERRA_API __declspec(dllimport)
#endif