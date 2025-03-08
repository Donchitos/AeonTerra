#pragma once

// Properly handle STATIC vs DLL exports
#ifdef AEONTERRA_CORE_STATIC
#define AEONTERRA_API
#else
#ifdef AEONTERRA_CORE_EXPORTS
#define AEONTERRA_API __declspec(dllexport)
#else 
#define AEONTERRA_API __declspec(dllimport)
#endif
#endif