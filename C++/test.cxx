// vector::operator[]
#include <iostream>
#include <vector>

void rd(std::vector<int> &df) {  
  for (std::vector<int>::iterator it = df.begin() ; it !=df.end(); ++it) {
    std::cout << *it << std::endl;
  }
}

int main ()
{
  std::vector<int> myvector (10);   // 10 zero-initialized elements
  std::vector<int>::size_type sz = myvector.size();

  // assign some values:
  for (unsigned i=0; i<sz; i++) myvector[i]=i;

  // reverse vector using operator[]:
  for (unsigned i=0; i<sz/2; i++)
  {
    int temp;
    temp = myvector[sz-1-i];
    myvector[sz-1-i]=myvector[i];
    myvector[i]=temp;
  }

  rd(myvector);

  std::cout << "myvector contains:";
  for (unsigned i=0; i<sz; i++)
    std::cout << ' ' << myvector[i];
  std::cout << '\n';

  return 0;
}
