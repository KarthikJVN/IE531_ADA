#include <iostream >
using name space std ;
void print integer (int num)
{
put char ( num % 10 + ’A ’ ) ;
if ( num / 10 )
print integer ( num / 10 ) ;
}
int main ( )
{
print integer ( 1234 ) ; cout << endl ;
}
