#ifndef STRLWR
#define STRLWR
#include <cctype>

inline char* lwr( char* str )
{
   char* orig = str;
   // process the string
   for ( ; *str != '\0'; str++ )
       *str = tolower(*str);
   return orig;
}
inline char* upr( char* str )
{
   char* orig = str;
   // process the string
   for ( ; *str != '\0'; str++ )
       *str = toupper(*str);
   return orig;
}
inline bool cmp( char* str1, char * str2 )
{
   // process the string
   for ( ; *str1 != '\0' && *str2 != '\0';)
   {
       if(*str1 != *str2) return 1;
       str1++;
       str2++;
   }
   return 0;
}
#endif
