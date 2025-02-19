// This is trying to construct the recurrence relation ri from ri-1
#include<iostream>
#include<vector>
#include<algorithm>
#include<cstdio>
using namespace std;
template <typename T>
void display(const std::vector<T>  &v);
template <typename T>
void display(const std::vector<T> &v){
    cout<<"{";
    for (size_t i = 0; i < v.size(); i++){
         cout << v[i] << ",";
    }
    cout<<"}"<<endl;
}
template <typename T>
std::vector<T> berlekamp_massey(const std::vector<T> &s){
    cout<<"In berlekamp_massey"<<endl;
    display(s);
    int n=(int)s.size();
    printf("n= %d \n",n);
    int l=0;// length of c- c.size()
    int m=0;
    vector<T> old_c(n),c(n);//old_c and c declared
    T delta=old_c[0]=c[0]=1;//delta, old_c[0] and c[0] initialized to 1.
    display(c);
    display(old_c);
    cout<<"delta= "<<delta<<endl;
    for (int i=0;i<n;i++){
        cout<<"_______________________________________________________________________"<<endl;
        cout<<"In for loop. i= "<<i<<endl;
        T diff=s[i];
        printf("diff= %d \n",diff);
        for (int j=1;j<=l;j++){// j=1 to c.size()
            cout<<"In middle for loop.LRR j= "<<j<<endl;
            diff+=c[j]*s[i-j];
        }
        cout<<"temp = "<<diff<<endl;
        if (diff==0){
            printf("NO failure\n");
            continue;
        
        }
        vector<T> temp_c=c;
        display(temp_c);
        T error_factor=diff/delta;
        cout<<"error_factor= "<<error_factor<<endl;
        for(int j=m;j<n;j++){
            cout<<"In last loop. j= "<<j<<endl;
            c[j]-=error_factor*old_c[j-m];
            cout<<"c["<<j<<"]= "<<c[j]<<endl;
        }
        if(2*l<=i){// DAFUQ is this? i=0 0=0, 
        cout<<"in DAFUQ is this"<<endl;
            l=i+1-l;// L is updated
            cout<<"updated l = "<<l<<endl;
            m=i;
            // m=0;
            old_c=temp_c;
            delta=diff;
        }
    }
    printf("=======================================================================\n");
    printf("End of for loop\n");
    cout<<"c before resizing"<<endl;
    display(c);
    c.resize(l+1);
    display(c);
    cout<<"c after resizing"<<endl;
    c.erase(c.begin());
    cout<<"c after erasing"<<endl;
    display(c);
    for (T &x:c)
        x=-x;
    return c;
}
int main(){
    vector<int> s={1,2,4,8,13,20,28,215,757,2186};
    vector<int> c=berlekamp_massey(s);
    for (int x:c)
        cout<<x<<" ";
    cout<<endl;
    return 0;
}