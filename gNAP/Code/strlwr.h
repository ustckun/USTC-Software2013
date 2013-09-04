#ifndef STRLWR
#define STRLWR
inline char* strlwr( char* str )
{
   char* orig = str;
   // process the string
   for ( ; *str != '\0'; str++ )
       *str = tolower(*str);
   return orig;
}
inline char* strupr( char* str )
{
   char* orig = str;
   // process the string
   for ( ; *str != '\0'; str++ )
       *str = toupper(*str);
   return orig;
}
inline bool strcmp( char* str1, char * str2 )
{
   // process the string
   for ( ; *str1 != '\0' && *str2 != '\0';){
       str1++;
       str2++;
       if(*str1 != *str2) return 0;
   }
   return 1;
}
#endif
