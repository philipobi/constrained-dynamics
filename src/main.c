#include <stdio.h>
#include <wchar.h>

int main()
{
    wchar_t a = L'𝒶';
    char *mystring = "hello";
    wchar_t out;
    for (char *p_char = mystring; *p_char != '\0'; p_char++) {
        out = a + (*p_char - 'a');
        wprintf(L"%lc", a);
    }
    return 0;
}