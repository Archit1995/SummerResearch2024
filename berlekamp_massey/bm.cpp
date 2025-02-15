// NOTE: use g++ to compile the code 
// Berlekamp Massey
#include <iostream>
#include<vector>
#include <algorithm>
#include<cstdio>
using namespace std;
// In cpp we need to declare functions before using them
int lrr(std::vector<int>c,int s[],int index);
std::vector<int> failure(std::vector<int>coeff,int* seq,int findex,int index, int del);
void display(std::vector<int>v);
void display(int* v,int size);
// function overloading

std::vector<int> add(std::vector<int> a,std::vector<int> b);
std::vector<int>pad(std::vector<int>coeff,int num_zeros);

std::vector<int>pad(std::vector<int>coeff,int num_zeros){
    cout<<"In pad"<<endl;
    for (int i = 0; i < num_zeros; i++)
        coeff.push_back(0);
    return coeff;
}
// c_new=c+temp_c
std::vector<int> add(std::vector<int> a,std::vector<int> b){
    std::vector<int> c_new;
    int size_diff=a.size()-b.size();
    if(size_diff>0){
        cout<<"Size of a > size of b, implies b needs padding "<<endl;
        b=pad(b,size_diff);
    }
    else if (size_diff<0){
        cout<<"Size of a < size of b implies a needs padding"<<endl;
        a=pad(a,-size_diff);
    }

        cout<<"Size of a = size of b"<<endl;
    for (int i = 0; i < a.size(); i++)
        c_new.push_back(a[i]+b[i]);
    return c_new;
}
void display(std::vector<int>v, char name[]){
    // cout<<"In display(std::vector<int>v) \n";
    for (int i = 0; i < v.size(); i++)
    {
        printf("%s[%d]= %d ",name,i,v[i]);
    }
    printf("\n");
}
void display(int* v,int size,char name[]){
    // printf("In display(int* v) \n");
    // cout<<"v.size()= "<<sizeof(v)<<endl;
    for (int i=0;i<size;i++){
        printf("%s[%d]= %d ",name,i,v[i]);
    }   
    printf("\n");
}

std::vector<int> failure(std::vector<int>coeff, int* seq, int findex,int index,int del){
    cout<<"In failure"<<endl;
    printf("findex= %d \n",findex);
    // Step 1 d=old_c. Note indexing for d starts from 1
    vector<int>d(coeff);
    // Step 2 d=-d
    std::transform(d.begin(),d.end(),d.begin(),[](int a){return a*-1;});
    // Step 3 d=[1,d], indexing starts at 1
    d.insert(d.begin()+1,1);// insert 1 at d[1] 
    printf("d.size() =%lu \n",d.size());

    for (int i = 0; i < d.size(); i++)
    {
        printf("d[%d]= %d \n",i,d[i] );
    }
    int t=lrr(d,seq,findex+1);
    cout<<"lrr(d,seq,f+1) = "<<t<<endl;
    // Step 4 del/t
    int error_factor=del/t;
    cout<<"error_factor= "<<error_factor<<endl;
    // Step 5 
    // std::transform(d.begin(),d.end(),d.begin(),[](int a){return a*error_factor;});
    for (int i=1;i<d.size();i++){
        d[i]=d[i]*error_factor;
    }
    // Step 6
    int offset=index-findex-1;
    if (offset >0){
        d.insert(d.begin()+1,offset,0);
        
    }

    
    //
    return d;

}

int lrr(std::vector<int> coeff,int* seq,int index){
    if(index==0){
        return 0;
    }
    else{
        printf("__________________________________ \n");
        cout<<"In lrr"<<endl;
        int sum=0;
        // cout<<"index= "<<index<<endl;
        cout<<"coeff.size()= "<<coeff.size()<<endl;
        for (int j = 1; j < coeff.size(); j++)
        {
            printf("coeff[%d]= %d \n",j,coeff[j]);  
            sum+=coeff[j]*seq[index-j];
            printf("sum = %d \n",sum);
        }
    return sum;
    }
}
int main(){
    cout<<"In main"<<endl;
    int N=10;
    int f=0;
int s[N]={1,2,4,8,13,20,28,215,757,2186};
display(s,N,"s");
vector<int>c;
vector<int>old_c;
vector<int>temp_c;
cout<<"old_c length= "<<old_c.size()<<endl;

// cout<<__cplusplus;
old_c.push_back(0);

// cout<<c.begin()<<endl;
// cout<<"c[0] "<<c[0]<<endl;
try{
    c.insert(c.begin(),0);
    cout<<"c[1] "<<c[1]<<endl;
    throw 20;
}
catch(int err){
    cout<<err<<endl;
}

for (int i = 0; i < N ; i++)
{

    // if(i==5){
    //     break;
    // }

    printf("_______________________________________________________________________\n");
    printf("In for loop. i= %d \n",i);
    // cout<<"s["<<i<<"]= "<<s[i]<<" "<<endl;
    // printf("old_c \n");
    display(old_c,"old_c");
    display(c,"c");
    int temp=lrr(c,s,i);
    cout<<"temp= "<<temp<<endl;
    int delta=s[i]-temp;
    cout<< "delta = "<<delta<<endl;
   
    if(delta!=0){
        printf("FAILURE at i= %d \n",i);
        if (i==0){
            c.insert(c.begin()+1,s[0]);
            f=i;
        }
        else{
        temp_c=failure(old_c,s,f,i,delta);
        f=i;
        cout<<"temp_c.size()= "<<temp_c.size()<<endl;
        display(temp_c,"temp_c");
        display(c,"c");
        old_c=c;
        c=add(c,temp_c);
        display(c,"c");

        }
    }
    else{
        continue;
    

    }
    
}

return 0;
}
// There exists 2 global sequences c and c_old and a local sequence d consrucuted from c_old
// c is the current sequence of coefficients c=c+d
// Note: we use c_old to construct d such that lrr(d,f+1)=delta, 
// Case 1: i=f 
// Case 2: i=f+1
// Case 3 : i>f+1



