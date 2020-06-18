# 复习目录

- 常用函数：STL等常用数据结构（树、栈、队列、优先队列、向量、散列表、集合、字符串）和函数

- 数学问题（最小公倍数、最大公因数、快速幂、进制转换、质因数、大数）

- #### 搜索（BFS & DFS）

- 图论（并查集、最小生成树、最短路径、拓扑排序、关键路径）

- #### 动态规划

# 常用函数

## memset-对数组中每个元素赋值0/-1

  ```c++
#include <string.h>
int a[5];
memset(a,0,sizeof(a));
  ```

## getchar/putchar && gets/puts

```c++
#include <stdio.h>
getchar(); //输入单个字符 有时候还可以用来吸收换行符
putchar(); //输出单个字符
gets(); //输入一行字符串，识别\n作为输入结束
puts(); //输出一行字符串，并输出\n 若是字符串末尾没有\0则会乱码
```

## string.h

```c++
#include <string.h>
strlen(str);
strcmp(str1,str2);
//比较字典序，小于返回负整数，等于返回0，大于返回正整数
strcpy(str1,str2); //把str2复制给str1
strcat(str1,str2); //把str2接到str1后面
```

## math.h

```c++
#include <math.h>
fabs(-12.56);  //12.56
floor(-5.2); //-6 向下（小）取整
ceil(-5.2); //5 向上（大）取整
double db = pow(2.0,3.0); //8.000000
sqrt(2.0); //1.414214
log(1.0); //0.000000 以自然对数为底的对数
//log(a)b=logb/loga  必须用换底公式
const double pi = acos(-1.0);
sin(pi*45/180);
cos(pi*45/180);
tan(pi*45/180);  //弧度制！！
asin(1);
acos(-1.0);
atan(0);
round(3.5); //4 四舍五入
```

## printf

```c++
#include <iostream>
#include <cstdio>
%5d; //宽为5的整数，超过5位按实际输出，不够5位右对齐输出
%05d; //宽为5的整数，超过5位按实际输出，不够5位前置补0右对齐
%5.2f; //宽为5的浮点数，小数点有2位，小数点占一位, 不够5位右对齐输出
%6.9s; //表示显示一个长度不小于6且不大于9的字符串。若大于9, 则第9个字符以后的内容将被删除。  
%ld;  //表示输出long整数    
%lf;  //表示输出double浮点数
%-7d;  //表示输出7位整数左对齐 空位补齐     
%-10s; //表示输出10个字符左对齐

/*按整型输出，默认右对齐*/
printf("%d\n",PrintVal);
/*按整型输出，补齐4位的宽度，补齐位为空格，默认右对齐*/
printf("%4d\n",PrintVal);
/*按整形输出，补齐4位的宽度，补齐位为0，默认右对齐*/
printf("%04d\n",PrintVal);

/*按16进制输出，默认右对齐*/
printf("%x\n",PrintVal);
/*按16进制输出，补齐4位的宽度，补齐位为空格，默认右对齐*/
printf("%4x\n",PrintVal);
/*按照16进制输出，补齐4位的宽度，补齐位为0，默认右对齐*/
printf("%04x\n",PrintVal);

/*按8进制输出，默认右对齐*/
printf("%o\n",PrintVal);
/*按8进制输出，补齐4位的宽度，补齐位为空格，默认右对齐*/
printf("%4o\n",PrintVal);
/*按照8进制输出，补齐4位的宽度，补齐位为0，默认右对齐*/
printf("%04o\n",PrintVal);
```

# 常用模板代码

## 日期转换

闰年判断

```c++
bool IsLeapYear(int year){
    return (year%4==0 && year%100!=0) || (year%400==0);
}
```

## 反序数

```c++
int Reverse(int x){
    int revx=0;
    while(x != 0){
        revx*=10;
        revx+=x%10;
        x/=10;
    }
    return revx;
}
```

## 最大公约数 & 最小公倍数

```c++
//最大公约数求法
int GCD(int a,int b){ //辗转相除法 a大b小
    if(b==0){
        return a;
    }else{
        return GCD(b,a%b);
    }
}
//最小公倍数求法
a * b / GCD(a,b);
```

## 质数和分解质因数

### 1. 质数判定：若到sqrt(n)的整数，所有正整数均不能整除n，则可以判断n为素数（小于2必定不是素数）

### 2. 某范围内的素数筛法

```c++
const int MAXN = 10001;
vector<int> prime; //保存质数
bool isPrime[MAXN]; //标记数组
void Initial(){
    //初始化
    for(int i=0;i<MAXN;++i){
        isPrime[i] = true;
    }
    isPrime[0]=false;
    isPrime[1]=false;
    //素数筛法
    for(int i = 2;i<MAXN;++i){
        if(!isPrime[i]){ //非质数跳过
            continue;
        }
        prime.push_back(i);
        for(int j=i*i;j<MAXN;j+=i){  //注意从i*i开始！！！
            isPrime[j]=false; //质数的倍数必为非质数
        }
    }
}
```

### 3. 分解质因数：先用素数筛法,然后不断试除

## 快速幂

1) 求指数的二进制数

2) 求底数幂然后累乘

3) 矩阵快速幂，初始化为单位矩阵

```c++
int FastExp(int a,int b,int mod){
    int answer = 1; //初始化为1
    while(b!=0){  //不断将b转换为二进制数
        if (b%2==1){ //若位为1，累乘a的2^k次幂
            answer *= a;
            //answer %= mod; //求后三位，依题目而定
        }
        b/=2; //b不断转换
        a*=a; //a不断平方
        //a%=mod; //求后三位，依题目而定
    }
    return a
```

## 同余模定理

```c++
(a*b)%c=((a%c)*(b%c))%c;
(a+b)%c=((a%c)+(b%c))%c;
```

## BFS

```c++
#include <iostream>
#include <queue>
#define SIZE 1000
using namespace std;
queue<int> q_x, q_y, q_pos;
int step[4][2] = {
    {1,0},{0,1},{-1,0},{0,-1}
};
int m, n, board[SIZE][SIZE], start_x, start_y, target_x, target_y;
bool vis[SIZE][SIZE];
void bfs() {
    bool flag = false;  //判断是否成功的标志
    while (!q_x.empty() && !q_y.empty()) {  //判断队是否为空
        for (int k = 0; k < 4; ++k) {  //先遍历
            int tx = q_x.front() + step[k][0];  //状态转移
            int ty = q_y.front() + step[k][1];
            if (tx < 1 || ty < 1 || tx > n || ty > m)  //判断边界
                continue;
            if (!board[tx][ty] && !vis[tx][ty]) { 
                vis[tx][ty] = true;  //标记
                q_x.push(tx);  //入队 
                q_y.push(ty);
                q_pos.push(q_pos.front() + 1);
            }
            if (tx == target_x&&ty == target_y) { //判断到达目标的条件
                flag = true;
                break;
            }
        }
        if (flag)
            break;
        q_x.pop();  //队首出队
        q_y.pop();
        q_pos.pop();
    }
    return;
}
int main()
{
    while (cin >> m >> n) {
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                cin >> board[i][j];
            }
        }
        cin >> start_x >> start_y >> target_x >> target_y;
        q_x.push(start_x);
        q_y.push(start_y);  //初始入队
        q_pos.push(0);
        vis[start_x][start_y] = true;
        bfs();
        cout << q_pos.back() << endl;
    }
    return 0;
}
```
## DFS

```c++
int check(参数)
{
    if(满足条件)
        return 1;
    return 0;
}
 
void dfs(int step)
{
        判断边界
        {
            相应操作
        }
        尝试每一种可能
        {
               满足check条件
               标记
               继续下一步dfs(step+1)
               恢复初始状态（回溯的时候要用到）
        }
}
```

## 并查集

1. 用于判断图是否为连通图
2. 求图的连通分量

```c++
#include <iostream>
#include <cstdio>
using namespace std;
const int MAXN = 1000;
int father[MAXN];  //父亲结点
int height[MAXN];  //结点高度  
void Initial(int n){  //初始化
    for(int i=0;i<=n;i++){
        father[i]=i;  //每个结点的父亲为自己
        height[i]=0;  //每个结点初始高度为0
    }
}
int Find(int x){  //查找根节点
    if(x!=father[x]){  //路径压缩
        father[x]=Find(father[x]);
    }
    return father[x];
}
void Union(int x,int y){  //合并集合
    x=Find(x);
    y=Find(y);
    if(x!=y){  //矮树作为高数的子树
        if(height[x]<height[y]){
            father[x]=y;
        }else if (height[y]<height[x]){
            father[y]=x;
        }else{
            father[y]=x;
            height[x]++;
        }
    }
    return ;
}
//先Initial 再 一个一个输入 Union
//Find(i)==i; 说明有一个集合，即连通分量
```

## 最小生成树

连通所有结点的树中权值最小的即为最小生成树，实现同样用了并查集的思想：

1. Initial 后按边的权值排序

2. 判断边的两个顶点是否同属一个集合 （用Find)

3. 不同集合则合并（Union），加上这条关键边的长度
```c++
#include <algorithm>
struct Edge{
    int from; //边的起点
    int to; //边的终点
    int length;  //边的长度
    friend bool operator < (const Edge& e1,const Edge& e2) const {
        return e1.length < e2.length;
    }
};
Edge edge[MAXN*MAXN]; //边集
/*上边并查集的所有内容*/
void Kruskal(int n,int edgeNum){ //克鲁斯卡尔算法，n为点数，edgeNum为边数
    Initial(n);
    sort(edge,edge+edgeNum); //按权值排序，重载运算符，无需重写cmp
    int sum=0; //初始化最小生成树长度
    for(int i=0;i<edgeNum;++i){  //遍历所有的边
        Edge current = edge[i];
        if (Find(current.from)!=Find(current.to)){ //边的两端结点不同属一个连通子集吗？
            Union(current.from,current.to);  //不同属的话就合并连通
            sum+=current.length;  //属于最小生成树的一部分了
        }
    }
    return sum;
}
```

## 大数

### 定义类

```c++
#include<iostream>
#include<cstring>
#include<cstdio>
#include<iomanip>
#include<algorithm>
using namespace std;

#define MAXN 9999
#define MAXSIZE 10
#define DLEN 4

class BigNum
{
    private:
        int a[1500];    //可以控制大数的位数
        int len;       //大数长度
    public:
        BigNum(){ len = 1;memset(a,0,sizeof(a)); }   //构造函数
        BigNum(const int);       //将一个int类型的变量转化为大数
        BigNum(const char*);     //将一个字符串类型的变量转化为大数
        BigNum(const BigNum &);  //拷贝构造函数
        BigNum &operator=(const BigNum &);   //重载赋值运算符，大数之间进行赋值运算

        friend istream& operator>>(istream&,  BigNum&);   //重载输入运算符
        friend ostream& operator<<(ostream&,  BigNum&);   //重载输出运算符

        BigNum operator+(const BigNum &) const;   //重载加法运算符，两个大数之间的相加运算
        BigNum operator-(const BigNum &) const;   //重载减法运算符，两个大数之间的相减运算
        BigNum operator*(const BigNum &) const;   //重载乘法运算符，两个大数之间的相乘运算
        BigNum operator/(const int   &) const;    //重载除法运算符，大数对一个整数进行相除运算

        BigNum operator^(const int  &) const;    //大数的n次方运算
        int    operator%(const int  &) const;    //大数对一个int类型的变量进行取模运算
        bool   operator>(const BigNum & T)const;   //大数和另一个大数的大小比较
        bool   operator>(const int & t)const;      //大数和一个int类型的变量的大小比较

        void print();       //输出大数
};
```
### 将一个int类型的变量转化为大数
```c++
BigNum::BigNum(const int b)     
{
	int c,d = b;
	len = 0;
	memset(a,0,sizeof(a));
	while(d > MAXN)
	{
		c = d - (d / (MAXN + 1)) * (MAXN + 1);
		d = d / (MAXN + 1);
		a[len++] = c;
	}
	a[len++] = d;
}
```
### 将一个字符串类型的变量转化为大数
```c++
BigNum::BigNum(const char*s)     
{
	int t,k,index,l,i;
	memset(a,0,sizeof(a));
	l=strlen(s);
	len=l/DLEN;
	if(l%DLEN)
		len++;
	index=0;
	for(i=l-1;i>=0;i-=DLEN)
	{
		t=0;
		k=i-DLEN+1;
		if(k<0)
			k=0;
		for(int j=k;j<=i;j++)
			t=t*10+s[j]-'0';
		a[index++]=t;
	}
}
```
### 拷贝构造函数
```c++
BigNum::BigNum(const BigNum & T) : len(T.len)  
{
	int i;
	memset(a,0,sizeof(a));
	for(i = 0 ; i < len ; i++)
		a[i] = T.a[i];
}
```
### 重载赋值运算符，大数之间进行赋值运算
```c++
BigNum & BigNum::operator=(const BigNum & n)   
{
	int i;
	len = n.len;
	memset(a,0,sizeof(a));
	for(i = 0 ; i < len ; i++)
		a[i] = n.a[i];
	return *this;
}
```
### 重载输入运算符
```c++
istream& operator>>(istream & in,  BigNum & b)   
{
	char ch[MAXSIZE*4];
	int i = -1;
	in>>ch;
	int l=strlen(ch);
	int count=0,sum=0;
	for(i=l-1;i>=0;)
	{
		sum = 0;
		int t=1;
		for(int j=0;j<4&&i>=0;j++,i--,t*=10)
		{
			sum+=(ch[i]-'0')*t;
		}
		b.a[count]=sum;
		count++;
	}
	b.len =count++;
	return in;

}
```
### 重载输出运算符
```c++
ostream& operator<<(ostream& out,  BigNum& b)   
{
	int i;
	cout << b.a[b.len - 1];
	for(i = b.len - 2 ; i >= 0 ; i--)
	{
		cout.width(DLEN);
		cout.fill('0');
		cout << b.a[i];
	}
	return out;
}
```
### 两个大数之间的相加运算
```c++
BigNum BigNum::operator+(const BigNum & T) const   
{
	BigNum t(*this);
	int i,big;      //位数
	big = T.len > len ? T.len : len;
	for(i = 0 ; i < big ; i++)
	{
		t.a[i] +=T.a[i];
		if(t.a[i] > MAXN)
		{
			t.a[i + 1]++;
			t.a[i] -=MAXN+1;
		}
	}
	if(t.a[big] != 0)
		t.len = big + 1;
	else
		t.len = big;
	return t;
}
```
### 两个大数之间的相减运算
```c++
BigNum BigNum::operator-(const BigNum & T) const   
{
	int i,j,big;
	bool flag;
	BigNum t1,t2;
	if(*this>T)
	{
		t1=*this;
		t2=T;
		flag=0;
	}
	else
	{
		t1=T;
		t2=*this;
		flag=1;
	}
	big=t1.len;
	for(i = 0 ; i < big ; i++)
	{
		if(t1.a[i] < t2.a[i])
		{
			j = i + 1;
			while(t1.a[j] == 0)
				j++;
			t1.a[j--]--;
			while(j > i)
				t1.a[j--] += MAXN;
			t1.a[i] += MAXN + 1 - t2.a[i];
		}
		else
			t1.a[i] -= t2.a[i];
	}
	t1.len = big;
	while(t1.a[len - 1] == 0 && t1.len > 1)
	{
		t1.len--;
		big--;
	}
	if(flag)
		t1.a[big-1]=0-t1.a[big-1];
	return t1;
}
```
### 两个大数之间的相乘运算
```c++
BigNum BigNum::operator*(const BigNum & T) const   
{
	BigNum ret;
	int i,j,up;
	int temp,temp1;
	for(i = 0 ; i < len ; i++)
	{
		up = 0;
		for(j = 0 ; j < T.len ; j++)
		{
			temp = a[i] * T.a[j] + ret.a[i + j] + up;
			if(temp > MAXN)
			{
				temp1 = temp - temp / (MAXN + 1) * (MAXN + 1);
				up = temp / (MAXN + 1);
				ret.a[i + j] = temp1;
			}
			else
			{
				up = 0;
				ret.a[i + j] = temp;
			}
		}
		if(up != 0)
			ret.a[i + j] = up;
	}
	ret.len = i + j;
	while(ret.a[ret.len - 1] == 0 && ret.len > 1)
		ret.len--;
	return ret;
}
```
### 大数对一个整数进行相除运算
```c++
BigNum BigNum::operator/(const int & b) const   
{
	BigNum ret;
	int i,down = 0;
	for(i = len - 1 ; i >= 0 ; i--)
	{
		ret.a[i] = (a[i] + down * (MAXN + 1)) / b;
		down = a[i] + down * (MAXN + 1) - ret.a[i] * b;
	}
	ret.len = len;
	while(ret.a[ret.len - 1] == 0 && ret.len > 1)
		ret.len--;
	return ret;
}
```
### 大数对一个int类型的变量进行取模运算
```c++
int BigNum::operator %(const int & b) const    
{
	int i,d=0;
	for (i = len-1; i>=0; i--)
	{
		d = ((d * (MAXN+1))% b + a[i])% b;
	}
	return d;
}
```
### 大数的n次方运算
```c++
BigNum BigNum::operator^(const int & n) const    
{
	BigNum t,ret(1);
	int i;
	if(n<0)
		exit(-1);
	if(n==0)
		return 1;
	if(n==1)
		return *this;
	int m=n;
	while(m>1)
	{
		t=*this;
		for( i=1;i<<1<=m;i<<=1)
		{
			t=t*t;
		}
		m-=i;
		ret=ret*t;
		if(m==1)
			ret=ret*(*this);
	}
	return ret;
}
```
### 大数和另一个大数的大小比较
```c++
bool BigNum::operator>(const BigNum & T) const   
{
	int ln;
	if(len > T.len)
		return true;
	else if(len == T.len)
	{
		ln = len - 1;
		while(a[ln] == T.a[ln] && ln >= 0)
			ln--;
		if(ln >= 0 && a[ln] > T.a[ln])
			return true;
		else
			return false;
	}
	else
		return false;
}
```
### 大数和一个int类型的变量的大小比较
```c++
bool BigNum::operator >(const int & t) const    
{
	BigNum b(t);
	return *this>b;
}
```
### 输出大数
```c++
void BigNum::print()    
{
	int i;
	cout << a[len - 1];
	for(i = len - 2 ; i >= 0 ; i--)
	{
		cout.width(DLEN);
		cout.fill('0');
		cout << a[i];
	}
	cout << endl;
}
```
### 使用
```c++
BigNum a,b,c;
char s1[1000];
char s2[1000];
int main(){
    scanf("%s",s1);
	scanf("%s",s2);

	a = s1;
	b = s2;
	c = a + b;
	cout<<a<<" + "<<b<<" = "<<c<<endl;      //重载了<<符号，只能用cout输出

	a = s1;
	b = s2;
    c = a - b;
    cout<<a<<" - "<<b<<" = "<<c<<endl;

    a = s1;
    b = s2;
    c = a * b;
    cout<<a<<" * "<<b<<" = "<<c<<endl;

    a = s1;
    c = a / 1000;
    cout<<a<<" / 1000 = "<<c<<endl;         //只能除int类型，得到的结果为大数类

    a = s1;
    b = a % 1000;                       //只能对int类型取模，得到的结果为大数类
    cout<<a<<" mod 1000 = "<<c<<endl;

    a = s1;
    c = a^10;
    cout<<a<<" ^ 10 = "<<c<<endl;           //大数的n次方，n在int范围内

    a = s1;
    b = s2;
    if(a>b) puts("yes");                //重载运算符，可以和大数类或int类型比大小
    if(a>1000)  puts("yeah");           //只能用>，不能用<、<=、.>=

    return 0;
}
```

## 最短路径
### Dijkstra
```c++
#include<iostream>
#include<cstdio>
#include<algorithm>
#include<utility>
#include<vector>
#include<queue>
#define MAXN 20000
#define INF 2147483647
using namespace std;
typedef pair<int,int> pii;
priority_queue<pii, vector<pii>, greater<pii> > pq;
struct edge
{
    int to;
    int cost;
};
vector<edge> G[MAXN];//g[i]--i to g[i].to cost cost
int n, m, s;
int dis[MAXN];
void dijk(int s)
{
    for(int i = 1; i <= n; i++)
        dis[i] = INF;
    dis[s] = 0;
    pq.push(make_pair(0,s));
   // cout<<dis[s]<<endl;
    while(!pq.empty())
    {
        pii u = pq.top();
        pq.pop();
        int x = u.second; // bian hao
        //cout<<x<<endl;
        for(int i = 0; i < G[x].size(); i++)
        {
            edge e = G[x][i];
            if(dis[e.to] > dis[x] + e.cost)
            {
                dis[e.to] = dis[x] + e.cost;
                pq.push(make_pair(dis[e.to], e.to));
               // cout<<dis[e.to]<<endl;
            }
        }
    }
}
int main()
{
    cin >> n >> m >> s;
    int from, to, cost;
    edge in;
    for(int i = 0; i < m; i++)
    {
        scanf("%d%d%d",&from ,&to ,&cost);
        in.to = to; in.cost = cost;
        G[from].push_back(in);
    }
   // cout<<endl;
    dijk(s);
    for(int i = 1; i <= n; i++)
        printf("%d ", dis[i]);
    return 0;
}
```
### Floyd
```c++
void Floyd(){
    for(int k=0;k<len;k++)
        for(int i=0;i<len;i++)
            for(int j=0;j<len;j++)
                map[i][j]=min(map[i][j],map[i][k]+map[k][j]);
}
```

## 关键路径

```c++
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
using namespace std;
struct node{
    int eEarlist;
    int eLatest;
    int eDiff;
    int expd;
    int startV, endV;
};
void DFS(int u);
int G[510][510];
const int MAX = 0x7FFFFFFF;
vector<int> vEarlist, vLatest;//状态最早开始时间和最迟开始时间
vector<node> actvt;//活动记录
vector<int> topoSeq;//拓扑排序的序列
vector<int> inDegree;//入度（配合拓扑排序使用）
vector<vector<int> > path;//关键路径（可能不止一条）
vector<int> tempPath;//DFS作临时路径
vector<vector<int> > pre;//状态的父状态（通过某条活动a->b，则pre[b]=a（其中一条））
int N, M, startV, endV;
int main()
{
    fill(G[0], G[0]+510*510, MAX);
    cin >> N >> M >> startV >>endV;
    vEarlist.resize(N); vLatest.resize(N);
    actvt.resize(M);
    inDegree.resize(N);
    fill(inDegree.begin(), inDegree.end(), 0);
    for(int i=0; i<M; i++){
        int u, v, expd;
        cin >> u >> v >> expd;
        G[u][v] = expd;
        actvt[i].expd = expd; actvt[i].startV = u; actvt[i].endV = v;
        inDegree[v]++;
    }
    //1.先进行拓扑排序，获得拓扑序列
    queue<int> q;
    q.push(startV);
    while(!q.empty()){
        int u = q.front();
        q.pop();
        for(int i=0; i<N; i++){
            if(G[u][i] != MAX){//存在这条路
                inDegree[i]--;
                if(inDegree[i]==0){
                    q.push(i);
                }
            }
        }
        topoSeq.push_back(u);
    }
    //2.依据拓扑序列，更新顶点状态：顶点最早开始时间
    fill(vEarlist.begin(), vEarlist.end(), 0);
    for(int i=0; i<topoSeq.size(); i++){
        int u = topoSeq[i];
        for(int j=0; j<N; j++){
            if(G[u][j] != MAX && vEarlist[u] + G[u][j] > vEarlist[j]){
                vEarlist[j] = vEarlist[u] + G[u][j];
            }   
        }
    }
    //3.找到最终结束点的最早开始时间，更新vLatest的最终点
    fill(vLatest.begin(), vLatest.end(), vEarlist[endV]);
    //4.对拓扑排序的逆序，进行操作，找到其相邻边，更新vLatest
    reverse(topoSeq.begin(), topoSeq.end());//直接改就行，因为用不到了
    for(int i=0; i<N; i++){
        int u = topoSeq[i];//此时为逆序，所以u在后
        for(int j=0; j<N; j++){
            if(G[j][u] != MAX && vLatest[u] - G[j][u] < vLatest[j]){//有到达u的点
                vLatest[j] = vLatest[u] - G[j][u];
            }
        }
    }
    //5.对所有活动，更新其最早开始时间-->顶点最早开始时间
    //              更新其最迟开始时间-->后顶点最迟开始时间-时长
    //              更新可休息时间-->最迟-最早
    //              如果可休息时间为0-->pre[后].push_back(前)
    //              
    pre.resize(M);
    for(int i=0; i<M; i++){
        actvt[i].eEarlist = vEarlist[actvt[i].startV];
        actvt[i].eLatest  = vLatest[actvt[i].endV]     -     actvt[i].expd;
        actvt[i].eDiff    = actvt[i].eLatest           -     actvt[i].eEarlist;
        if(actvt[i].eDiff == 0){
            pre[actvt[i].endV].push_back(actvt[i].startV);
        }
    }     
    //6.对结尾进行DFS，获得关键路径
    // DFS(endV);
    // path存储路径
    DFS(endV);
    for(int i=0; i<path.size(); i++){
        for(int j=path[i].size()-1; j>=0; j--){
            if(j!=path[i].size()-1)cout << "-->";
            cout << path[i][j];
        }
        cout << endl;
    }
}
void DFS(int u)
{
    tempPath.push_back(u);
    if(u == startV){
        path.push_back(tempPath);
        tempPath.pop_back();
        return;
    }
    for(int i=0; i<pre[u].size(); i++){
        DFS(pre[u][i]);
    }
    tempPath.pop_back();
}
```

## 动态规划
### 求最长递增子序列的长度 O(N^2)
```c++
//求最长递增子序列的长度O(N^2)
int Arr[30010],List[30010];
int LIS(int *Arr,int N)    //arr[]存放的是待求数组
{
    int Max = 0;        //max为最大递增子序列的长度
    for(int i = 1; i <= N; ++i)
        List[i] = 1;    //lis[i] 存放i之前的最大递增子序列的长度，初始都为1

    for(int i = 2; i <= N; ++i)
        for(int j = 1; j < i; ++j)    //遍历i之前所有位置
            if(Arr[i] >= Arr[j] && List[i]<List[j]+1)
                List[i] = List[j] + 1;
            //arr[i]>arr[j]为递增
            //lis[i]<lis[j] + 1确保为当前最长递增序列

    for(int i = 1; i <= N; ++i)
        if(Max < List[i])
            Max = List[i];

    return Max;
}
```
### 求最长递增子序列的长度 O(NlogN)
```c++
//求最长递增子序列的长度O(NlogN)
int Arr[10010],List[10010];
int Stack[10010];
int LIS(int *Arr,int N)
{
    int top = 0;
    Stack[top] = -1;
    for(int i = 1; i <= N; ++i)
    {
        if(Arr[i] > Stack[top])
            Stack[++top] = Arr[i];
        else
        {
            int low = 1;
            int high = top;
            while(low <= high)
            {
                int mid = (low + high)/2;
                if(Arr[i] > Stack[mid])
                    low = mid + 1;
                else
                    high = mid - 1;
            }
            Stack[low] = Arr[i];
        }
    }
    return top;
}
```
### 最大连续子序列和
给定一个整数序列，找出所有连续子序列中元素和最大的一个，并找到起点和终点。
```c++
int a[110000],N,pos1,pos2,Start,End;
//Start、End存储最大连续子序列的起点和终点
int MaxSubSum(int *a)
{
    int MaxSum = a[0],Sum = a[0];
    pos1 = pos2 = Start = End = 0;
    for(int i = 1; i < N; ++i)
    {
        Sum += a[i];
        if(Sum < a[i])
        {
            Sum = a[i];
            pos1 = i;
            pos2 = i;
        }
        else
        {
            pos2 = i;
        }

        if(MaxSum < Sum)
        {
            MaxSum = Sum;
            Start = pos1;
            End = pos2;
        }
    }
    return MaxSum;
}
```
### 最大连续子矩阵和
给你一个N，接下来是N*N的矩阵。数有正有负，求最大的子矩阵和。
```c++
#include<stdio.h>
#include<string.h>
#include<algorithm>
using namespace std;

int map[110][110],dp[110][110];
int main()
{
    int N,a;
    while(~scanf("%d",&N) && N)
    {
        memset(map,0,sizeof(map));
        memset(dp,0,sizeof(dp));
        for(int i = 1; i <= N; i++)
        {
            for(int j = 1; j <= N; j++)
            {
                scanf("%d",&a);
                map[i][j] = map[i][j-1] + a;
                //map[i][j]表示第i行前j列的和
            }
        }
        int Max = -0xffffff0;
        for(int j = 1; j <= N; j++)
        {
            for(int i = 1; i <= j; i++)
            {
                dp[i][j] = 0;
                for(int k = 1; k <= N; k++)
                {
                    //ans求的是前k行，第i到第j列的最大和
                    dp[i][j]= max(dp[i][j]+map[k][j]-map[k][i-1],map[k][j]-map[k][i-1]);
                    if(dp[i][j] > Max)
                        Max = dp[i][j];
                }
            }
        }
        printf("%d\n",Max);
    }
    return 0;
}
```
### 最大M个连续子段的和
```c++
#include<stdio.h>
#include<string.h>
#include<algorithm>
using namespace std;

int dp[1000010];
int maxn[1000010];
int num[1000010];
int main()
{
    int M,N;
    while(~scanf("%d%d",&M,&N))
    {
        dp[0] = maxn[0] = 0;
        for(int i = 1; i <= N; i++)
        {
            scanf("%d",&num[i]);
            dp[i] = maxn[i] = 0;
        }
        int MAXN;
        for(int i = 1; i <= M; i++)//分为i段
        {
            MAXN = -0xffffff0;
            for(int j = i; j <= N; j++)//第j个数字
            {
                dp[j] = max(dp[j-1]+num[j],maxn[j-1]+num[j]);
                maxn[j-1] = MAXN;
                MAXN = max(MAXN,dp[j]);
            }
        }
        printf("%d\n",MAXN);
    }

    return 0;
}
```
### 最大不连续子序列和
给你一个矩阵，不能选择每行中相邻的数字，也不能选当前行的上一行和下一行，问使所选数和最大的值是多少？
对于每一行，都是求最大不连续子段和。
```c++
#include<stdio.h>
#include<string.h>
#include<algorithm>
using namespace std;
const int MAXN = 200000;
int dpa[MAXN+20],dpb[MAXN+20],row[MAXN+20];
int main()
{
    int M,N,num;
    while(~scanf("%d%d",&M,&N))
    {
        memset(row,0,sizeof(row));
        for(int i = 0; i < M; i++)
        {
            dpa[0] = dpb[0] = 0;
            for(int j = 0; j < N; j++)
            {
                scanf("%d",&num);
                dpa[j+1] = max(dpa[j],dpb[j]);// dp[j+1] 是到j为止，不吃j所能吃到的最大值
                dpb[j+1] = dpa[j] + num;//吃j所能吃到的最大值
            }
            row[i] = max(dpa[N],dpb[N]);
        }
        dpa[0] = dpb[0] = 0;
        for(int i = 0; i < M; i++)
        {
            dpa[i+1] = max(dpa[i],dpb[i]);
            dpb[i+1] = dpa[i] + row[i];
        }
        int ans = max(dpa[M],dpb[M]);
        printf("%d\n",ans);
    }
    return 0;
}
```
### 最长公共子序列
给定两个序列，找出在两个序列中同时出现的最长子序列的长度。一个子序列是出现在相对顺序的序列，但不一定是连续的。
#### 1.求最长公共子序列长度
```c++
char s1[220],s2[220];
int dp[220][220];
//求串s1和串s2的公共子序列
int lcs(char *s1,char *s2)
{
    int len1 = strlen(s1);
    int len2 = strlen(s2);
    for(int i = 0; i <= len1; ++i)
    {
        for(int j = 0; j <= len2; ++j)
        {
            if(i == 0 || j == 0)
                dp[i][j] = 0;
            else if(s1[i-1] == s2[j-1])
                dp[i][j] = dp[i-1][j-1] + 1;
            else
                dp[i][j] = max(dp[i-1][j],dp[i][j-1]);
        }
    }
    return dp[len1][len2];
}
```
#### 2.求最长公共子序列长度，并输出路径
```c++
int dp[110][110],pre[110][110],len1,len2;
char s1[110],s2[110];

void LCS(char *s1,char *s2)
{
    for(int i = 0; i <= len1; ++i)
        pre[i][0] = 1;
    for(int i = 0; i <= len2; ++i)
        pre[0][i] = 2;
    //得到最长公共子序列，并标记dp[i][j]的上一个状态，用来回溯寻找路径
    for(int i = 0; i <= len1; ++i)
    {
        for(int j = 0; j <= len2; ++j)
        {
            if(i == 0 || j == 0)
                dp[i][j] = 0;
            if(s1[i-1] == s2[j-1])
            {
                dp[i][j] = dp[i-1][j-1] + 1;
                pre[i][j] = 0;
            }
            else if(dp[i-1][j] >= dp[i][j-1])
            {
                dp[i][j] = dp[i-1][j];
                pre[i][j] = 1;
            }
            else
            {
                dp[i][j] = dp[i][j-1];
                pre[i][j] = 2;
            }
        }
    }
}

void Print(int i,int j) //回溯输出新的字符串序列
{
    if(i == 0 && j == 0)
        return ;
    if(pre[i][j] == 0)
    {
        Print(i-1,j-1);
        printf("%c",s1[i-1]);
    }
    else if(pre[i][j] == 1)
    {
        Print(i-1,j);
        printf("%c",s1[i-1]);
    }
    else
    {
        Print(i,j-1);
        printf("%c",s2[j-1]);
    }
}

void solve(char *s1,char *s2)
{
    len1 = strlen(s1);
    len2 = strlen(s2);
    LCS(s1,s2);
    Print(len1,len2);
    printf("\n");
}
```
### 最长回文子序列
给一个字符串，找出它的最长的回文子序列LPS的长度。例如，如果给定的序列是“BBABCBCAB”，则输出应该是7，“BABCBAB”是在它的最长回文子序列。
```c++
char s[2020];
int dp[2020][2020];
//dp[i][j]表示s[i~j]最长回文子序列
int LPS(char *s)
{
    memset(dp,0,sizeof(dp));
    int len = strlen(s);
    for(int i = len-1; i >= 0; --i)
    {
        dp[i][i] = 1;
        for(int j = i+1; j < len; ++j)
        {
            if(s[i] == s[j])    //头尾相同，最长回文子序列为去头尾的部分LPS加上头和尾
                dp[i][j] = dp[i+1][j-1] + 2;
            else    //头尾不同，最长回文子序列是去头部分的LPS和去尾部分LPS较长的
                dp[i][j] = max(dp[i][j-1],dp[i+1][j]);
        }
    }

    return dp[0][len-1];
}
```
### 最长回文子串
给一个字符串，找出它的最长的回文子串(连续的串)的长度。
```c++
char str[2000020],s[2000020];
//str为待求的原串,s为处理后的新串
int P[2000020];
//P［i］记录的是以字符str［i］为中心的最长回文串的半径
void Pre(char *str)
{
    int len = strlen(str);
    s[0] = '$';
    s[1] = '#';
    for(int i = 0; i < len; ++i)
    {
        s[i*2+2] = str[i];
        s[i*2+3] = '#';
    }
    s[len*2+2] = '\0';
}
//返回最长回文子串长度
int Manacher(char *s)
{
    int Max = 0;
    int len = strlen(s);
    int id = 0;
    for(int i = 1; i < len; ++i)
    {
        if(Max > i)
            P[i] = min(P[2*id-i],P[id]+id-i);
        else
            P[i] = 1;
        while(s[i+P[i]] == s[i-P[i]])
            P[i]++;
        if(P[i]+i > Max)
        {
            Max = P[i]+i;
            id = i;
        }
    }
    int ans = 0;
    for(int i = 2; i < len; ++i)
        if(ans < P[i])
            ans = P[i];
    return ans-1;//返回最长回文子串长度
}


int main()
{
    while(~scanf("%s",str))
    {
        Pre(str);
        printf("%d\n",Manacher(s));
    }

    return 0;
}
```
### 最小编辑距离
给定一个长度为m和n的两个字符串，设有以下几种操作：替换（R），插入（I）和删除（D）且都是相同的操作。寻找到转换一个字符串插入到另一个需要修改的最小（操作）数量。
```c++
int dp[1010][1010],len1,len2;
char s1[1010],s2[1010];
int EditDist(char *s1,char *s2)
{
    int len1 = strlen(s1);
    int len2 = strlen(s2);
//当两个字符串的大小为0，其操作距离为0。
//当其中一个字符串的长度是零，需要的操作距离就是另一个字符串的长度. 
    for(int i = 0; i <= len1; ++i)
        dp[i][0] = i;
    for(int i = 0; i <= len2; ++i)
        dp[0][i] = i;

    for(int i = 1; i <= len1; ++i)
    {
        for(int j = 1; j <= len2; ++j)
        {
            if(s1[i-1] == s2[j-1])  //对齐s1[i-1]和s2[j-1]，不需改变
                dp[i][j] = dp[i-1][j-1];
            else    
                dp[i][j] = min(dp[i-1][j],min(dp[i][j-1],dp[i-1][j-1])) + 1;
            //s1前缀右对齐,s2前缀右为' ',删除s1第i个字符 -> dp[i-1][j]
            //s2前缀右对齐,s1前缀右为' ',删除s2第j个字符 -> dp[i][j-1]
            //s1前缀右对齐,s2前缀右对齐,i、j不一样，替换 -> dp[i-1][j-1]
        }
    }
    return dp[len1][len2];
}
```
### 01背包
有 N 件物品和一个容量为 V 的背包。放入第 i 件物品耗费的空间是 C i ，得到的价值是 W i 。求解将哪些物品装入背包可使价值总和最大。
要求恰好装满背包，在初始化时 F[0] = 0 ，其它 F[1~V ] 均设为 −∞ 。
没有要求必须把背包装满，而是只希望价格尽量大，初始化时应该将 F[0~V ] 全部设为 0 。
求方案数时，在初始化时F[0] = 1，其他F[1~V]均设为0。

```c++
int c[1100],w[1100],dp[1100],V;
//c[]：物品所占容量;w[]物品的价值;V为背包容量
memset(dp,0,sizeof(dp));
for(int i = 0; i < N; ++i)//第i件物品
{
     for(int j = V; j >= c[i]; j--)//填满空间j
     {
         dp[j] = max(dp[j],dp[j-c[i]] + w[i]);
     }
}
```
### 完全背包
有 N 种物品和一个容量为 V 的背包，每种物品都有无限件可用。放入第 i 种物品的耗费的空间是 C i ，得到的价值是 W i 。求解：将哪些物品装入背包，可使这些物品的耗费的空间总和不超过背包容量，且价值总和最大。
```c++
int c[1100],w[1100],dp[1100],V;
//c[]：物品所占容量;w[]物品的价值;V为背包容量
memset(dp,0,sizeof(dp));
for(int i = 0; i < N; ++i)
{
    for(int j = c[i]; j <= V; j++)
    {
        dp[j] = max(dp[j],dp[j-c[i]] + w[i]);
    }
}
```
### 多重背包
有 N 种物品和一个容量为 V 的背包。第 i 种物品最多有 M i 件可用，每件耗费的空间是 C i ，价值是 W i 。求解将哪些物品装入背包可使这些物品的耗费的空间总和不超过背包容量，且价值总和最大。
二进制优化做法：
```c++
int c[110],w[110],m[110],dp[100010],V;
//c[]:物品所占容量;w[]物品的价值;m[]物品的数量;V为背包容量
void ZeroOne(int cost,int weight)//01背包
{
    for(int i = V; i >= cost; i--)
        dp[i] = max(dp[i],dp[i-cost]+weight);
}

void Complete(int cost,int weight)//完全背包
{
    for(int i = cost; i <= V; i++)
        dp[i] = max(dp[i],dp[i-cost]+weight);
}

void Multiple(int cost,int weight,int cnt)//多重背包
{
//如果总容量比这个物品的容量要小，那么这个物品可以直接取完，相当于完全背包
    if(V <= cnt*cost)
    {
        Complete(cost,weight);
        return;
    }
    else//否则就将多重背包转化为01背包
    {
        int k = 1;
        while(k <= cnt)
        {
            ZeroOne(k*cost,k*weight);
            cnt -= k;
            k <<= 1;
        }
        ZeroOne(cnt*cost,cnt*weight);
    }
}
    /*
    for(int i = 0; i <= V; i++)//初始化,不要求恰好装满背包
        dp[i] = 0;  
    */
    for(int i = 0; i <= V; i++)//初始化：是否恰好装满背包
        dp[i] = -0xffffff0;
    dp[0] = 0;
    for(int i = 0; i < N; i++)
        Multiple(v[i],v[i],c[i]);
```
### 二维费用背包
```c++
int c1[110],c2[110],w[110],dp[1010][1010],V1,V2;
memset(dp,0,sizeof(dp));  

for(int i = 0; i < N; i++)//第i个  
{  
    for(int j = c1[i]; j <= V1; j++)//一维费用 
    {  
        for(int k = c2[i]; k <= V2; k++)//二维费用  
        {  
            dp[j][k] = max(dp[j][k],dp[j-c1[i]][k-c2[i]] + w[i]);  
        }  
    }  
}  
```
### 切割钢条
给定一段长度为n英寸的钢条和一个价格表Pi，求切割方案，使得销售收益Rn最大。
自底向上法  
```c++
#include<stdio.h>  
#include<algorithm>  
using namespace std;  
const int INF = 0xffffff0;  
int p[110],r[110];//r[n]来保存子问题  

int BOTTOM_UP_CUT_ROD(int n)  
{  
    r[0] = 0;//长度为0的钢条没有收益  
    for(int j = 1; j <= n; j++)//对j=1,2,3,…,n按升序求解每个规模为j的子问题。  
    {  
        int q = -INF;  
        for(int i = 1; i <= j; i++)  
        {  
            q = max(q,p[i]+r[j-i]);//直接访问数组r[j-i]来获得规模为j-i的子问题的解  
        }  
        r[j] = q;  
    }  
    return r[n];  
}  
int main()  
{  
    int N;  
    while(~scanf("%d",&N))  
    {  
        for(int i = 1; i <= N; i++)  
            scanf("%d",&p[i]);  

        int ans = BOTTOM_UP_CUT_ROD(N);  
        printf("%d\n",ans);  
    }  
    return 0;  
}  
```
### 最大矩形问题
给你一个直方图，告诉你各个条形矩形的高度，求基线对齐构成的矩形中面积最大的矩形的面积。
```c++
#include<stdio.h>
#include<string.h>

int l[100010],r[100010];
__int64 h[100010];
int main()
{
    int N;
    while(~scanf("%d",&N) && N!=0)
    {
        memset(h,0,sizeof(h));
        for(int i = 1; i <= N; i++)
        {
            scanf("%I64d",&h[i]);
            l[i] = r[i] = i;
        }

        l[0] = 1;
        r[N+1] = N;
        h[0] = -1;
        h[N+1] = -1;
        //这上边不加就会超时，不加的话下边就可能一直while，跳不出循环
        for(int i = 1; i <= N; i++)
        {
            while(h[l[i]-1] >= h[i])//找位置i的左边界
                l[i] = l[l[i]-1];
        }
        for(int i = N; i >= 1; i--)
        {
            while(h[r[i]+1] >= h[i])//找位置i的右边界
                r[i] = r[r[i]+1];
        }
        __int64 MaxArea = -0xffffff0;
        for(int i = 1; i <= N; i++)
        {
            if(h[i]*(r[i]-l[i]+1) > MaxArea)
                MaxArea = h[i]*(r[i]-l[i]+1);
        }
        printf("%I64d\n",MaxArea);
    }
    return 0;
}
```

## 二叉树

*二叉排序树（左小右大）中序遍历为升序序列

### 定义&创建

```c++
struct Node{
    char data;
    Node * left;
    Node * right;
    Node(char c): data(c),left(NULL),right(NULL){} //构造函数写法
};
//创建二叉树
Node * root = Build(str);
Node * Build(string str){
    //判断返回空树条件
    if(str.size()==0){
        return NULL;
    }
    //创建新节点
    Node * root = new Node();
    //递归创建左右子树
    root -> left = Build(str);
    root -> right = Build(str);
    return root;
}
//插入节点也可以采用递归的方式，同样需要给出退出递归的判断条件，通常是判断为空（二叉排序树的建法），插入的话边输入边插
Node* Insert(Node* root,int x){
    if (root==NULL){
        root = new Node(x);
    }else if(x<root->data){
        root->left = Insert(root->left,x);
    }else{
        root->right=Insert(root->right,x);
    }
    return root;
}
```

### 遍历

  ```c++
//前序遍历
void PreOrder(Node* root){
    if (root == NULL){
        return;
    }
    Visit(root->data); //可以是任何的操作，printf之类的
    PreOrder(root->left);
    PreOrder(root->right);
    return;
}
//中序遍历
void InOrder(Node* root){
    if (root == NULL){
        return;
    }
    InOrder(root->left);
      Visit(root->data); //可以是任何的操作，printf之类的
      InOrder(root->right);
      return;
}
//后序遍历
void PostOrder(Node* root){
    if (root == NULL){
        return;
    }
    PostOrder(root->left);
    PostOrder(root->right);
    Visit(root->data); //可以是任何的操作，printf之类的
    return;
}
  ```

  # STL 总结

  ##  vector（变长数组）

  ```c++
#include <vector>
using namespace std
vector<int> vi;
vi.push_back(1);
vi.pop_back();
vi.size();
vi.clear();
vi.begin(); //vi[0]
vi.end(); //尾元素地址的下一个地址
vi.insert(vi.begin()+2,100);
vi.erase(vi.begin()+3);
vi.erase(vi.begin()+1,vi.begin()+3);  //[first,last)
for(vector<int>::iterator it = vi.begin(); it!= vi.end(); it++){
    printf("%d ",*it);  //*(it + 2)
}
  ```

  ## set

  ```C++
#include <set>    
using namespace std;
set<int> si;   //自动去重并按升序排序
si.insert(3);
  
set<int>::iterator it = si.find(2);
printf("%d\n",*it);  //*(si.find(2));

si.erase(si.find(100));
si.erase(100);
set<int>::iterator it = si.find(30);
si.erase(it,si.end());
  
si.size();
si.clear();
  ```

  ## string	

  ```c++
#include <string>
using namespace std;
string str;
string str = "abcd";
str.length();
str.size();
cin>>str;
cout<<str;
printf("%s\n",str.c_str());
for(string::iterator it = str.begin();it != str.end();it++){
    printf("%c",*it);
}
str = str1 + str2;
str1 < str2； //字典序比较
str.insert(3,str2); // 在str[3]处插入,即str[3]为插入起点
str.insert(str.begin()+3,str2.begin(),str2.end()); //在str[3]处插入str2 [begin,end)之间的字符串
str.erase(str.begin()+4);
str.erase(str.begin()+2,str.end()-1); 
str.erase(3,2); //3号位置开始的2个字符
str.clear();
str.substr(14,4); //14号位开始的4长度字符
if (str.find(str2) != string::npos){  //返回第一次出现的位置
    cout<<str.find(str2)<<endl;
    cout<<str.find(str2,7)<<endl; //从7号位开始匹配str2
}
str.replace(10,4,str2); //10号位开始，长度为4的字串替换为str2
str.replace(str.begin(),str.begin()+5,str3);// [0,5)范围字串替换为str3
  ```

  ## map

  ```c++
#include <map>
using namespace std;
map<string, int> mp;  // <key,value>
for(map<char, int>::iterator it = mp.begin();it !=mp.end();it++){
    printf("%c %d\n",it -> first, it -> second);
}  //自动实现从小到大的字典排序功能
map<char, int>::iterator it = mp.find('b');
mp.erase(it);
mp.erase('b');
mp.erase(it,mp.end()); //删除[it,end)之间的映射
mp.size();
mp.clear();
  ```

  ## queue

  ```c++
#include <queue>
using namespace std;
queue<int> q;
q.push(i);
q.front();
q.back();
q.pop(); //队首元素出队
q.empty() == true;//判断是否为空
q.size();
  ```

  ## priority_queue

  ```c++
#include <queue>
using namespace std;
priority_queue<int> q; //队首（堆顶）返回优先级最大的
priority_queue< int, vector<int>, less<int> > q; //与上面等价，数字越大优先级越大
priority_queue< int, vector<int>, greater<int> > q; //数字小的优先级大
q.top();  // 使用之前先用empty()判断优先队列是否为空
q.push(3);
q.pop();
q.empty();
q.size();
struct fruit{
    string name;
    int price;
    friend bool operator < (fruit f1, fruit f2){
        return f1.price > f2.price;
    }
}; //水果价格低的为优先级高
//或者用以下这种形式
struct fruit{
    string name;
    int price;
}f1,f2,f3;
struct cmp{
    bool operator()(fruit f1,fruit f2){
        return f1.price > f2.price; //与上面一致
    }
};
priority_queue<fruit, vector<fruit>, cmp> q;
//对于字符串或者数组比较重载使用 const 和 &
friend bool operator < (const fruit &f1, const fruit &f2){
    return f1.price > f2.price;
}
bool operator () (const fruit &f1,const fruit &f2){
    return f1.price >f2.price;
}
  ```

  ## stack

  ```c++
#include <stack>
using namespace std;
stack<int> st;
st.push(1);
st.top();
st.pop();
st.empty();
st.size();
  ```

  ## pair

  相当于有两个元素的结构体，用于代替二元结构体及其构造函数

  ```c++
#include <utility>  //#include <map> 也行
using namespace std;
pair<string, int> p;
pair<string, int> p("haha",5);
//临时创建pair
pair<string, int> p;
p = make_pair("haha",5);
p = pair<string, int> ("haha",5);
p.first = "haha";
p.second = 5;
//比较时候先以first大小为准，相等时比较second大小
p1 < p2;
//用途
map<string, int> mp;
mp.insert(make_pair("heihei",5));
mp.insert(pair<string, int>("haha",10));
  ```

  ## algorithm

  ```c++
#include <algorithm>
using namespace std;
  ```

  ```c++
max(x,y);
  ```

  ```c++
min(x,y);
  ```

  ```
abs(x);
  ```

  ```c++
int x=1,y=2;
swap(x,y);
  ```

  ```c++
int a[5]={0,1,2,3,4};
reverse(a,a+3);  //[0,3)
string str="abcdefghi";
reverse(str.begin()+2,str.begin()+6); //[2,6)
  ```

  ```c++
int a[10]={1,2,3};
do{
    printf("%d%d%d\n",a[0],a[1],a[2]);
}while(next_permutation(a,a+3)); 
/*给出一个序列a在全排列中的下一个序列，(a,a+3)就是定义这个序列a
    123
    132
    213
    ...
    321
*/
  ```

  ```c++
int a[5]={1,2,3,4,5};
fill(a,a+5,233); //[0,5)赋值为233
  ```

  ```c++
bool cmp(int a, int b){
    return a > b; //降序输出
}
sort(a,a+4,cmp); // cmp是比较函数，默认递增，a是数组
//vector、string、deque才可以用sort，set、map本身有序不行
sort(vi.begin(),vi.end(),cmp)
  ```

  ```c++
int a[10]={1,2,2,3,3,3,4,4,4,5,5,5};
//使用lower_bound()或upper_bound()要求用在有序的数组或容器
printf("%d, %d\n",lower_bound(a,a+10,3)-a); //-a之后可以直接获取位置索引
//用来寻找[0,10)范围内第一个值大于或等于3的元素的位置，返回指针或迭代器
printf("%d, %d\n",upper_bound(a,a+10,3)-a);//-a之后可以直接获取位置索引
//用来寻找[0,10)范围内第一个值大于3的元素的位置，返回指针或迭代器
//找不到则返回可以插入该元素的位置指针或迭代器，即元素应在的位置
  ```

## STL 对比

|          | set  | string | map  | vector    | queue | priority_queue | stack |
| -------- | ---- | :----: | ---- | --------- | ----- | -------------- | ----- |
| 尾部添加 |      |        |      | push_back | push  | push           | push  |
| pop      |      |        |      | pop_back  | pop   | pop            | pop   |
| size     | √    |   √    | √    | √         | √     | √              | √     |
| clear    | √    |   √    | √    | √         |       |                |       |
| begin    | √    |   √    | √    | √         | front | top            | top   |
| end      | √    |   √    | √    | √         | back  |                |       |
| insert   | √    |   √    | √    | √         |       |                |       |
| erase    | √    |   √    | √    | √         |       |                |       |
| empty    | √    |   √    | √    | √         | √     | √              | √     |
| find     | √    |   √    | √    |           |       |                |       |

# 参考资料
[1] 胡凡，曾磊主编.算法笔记[M].北京：机械工业出版社.2016.

[2] 王道论坛组编.计算机考研 机试指南[M].北京：电子工业出版社.2019.

[3] https://www.cnblogs.com/ray-coding-in-rays/p/6498345.html

[4] https://blog.csdn.net/u011676797/java/article/details/45424703

[5] https://blog.csdn.net/mig_davidli/article/details/8591504

[6] https://www.jianshu.com/p/a85fbd6527a7

[7] https://www.cnblogs.com/bestsort/p/10588907.html

[8] https://blog.csdn.net/Mashiro_ylb/java/article/details/78287724

[9] https://blog.csdn.net/qq_33904395/java/article/details/103794896
