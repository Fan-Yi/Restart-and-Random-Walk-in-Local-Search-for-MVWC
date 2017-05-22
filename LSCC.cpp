// =====================================================================================
//
//       Filename:  LSCC+BMP.cpp
//
//    Description:  This is a solver for weighted maximum clique problem based on SCC and BMP
//
//        Version:  1.0
//        Created:  2016.01
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Yiyuan Wang          yiyuanwangjlu@126.com
//                  Shaowei Cai          caisw@ios.ac.cn
//                  Minghao Yin          ymh@nenu.edu.cn

//         For example：./LSCC+BMP bio-yeast.mtx 1000 1 
// =====================================================================================

      
      
#include<limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include<sys/times.h>
#include<unistd.h>
#include <time.h>
#include <ctime>
//#include <vector>
#include <string>
#include<math.h>
#include <climits>
#include<assert.h>

#include "cliqueHash.h"
#include "myBijection.h"

using namespace std;

#define clique_hash_mode
//#define local_optimum_revisit_restart_mode
#define first_improved_revisit_restart_mode

//#define visited_vertex_set_mode

#define rand_drop_mode

//#define dynamic_tabu_management_mode
//#define aspiration_mode

//#define promising_starting_neighborhood_mode

//#define strategy_analysis_mode

//#define debug_mode


int BMP = INT_MAX;
//int BMP = 100;

tms start, finish;
double time_limit;
double real;
int lbest;
     //int M_iter = 0;
struct Edge1{
	int v1;
	int v2;
};

Edge1 *edge;  //MAXE
int *v_degree_tmp;//MAXV

int *adjaclen;//MAXV

double real_solve1=-1;
double real_solve2=-1;
int **neighbor;
int *neighbor_len;//MAXV
int *conf_change;//MAXV

int *time_stamp;//MAXV




int *temp_array;//MAXV

char * File_Name;

int *vectex;//MAXV
int *funch;//MAXV
int *address;//MAXV
int *tabuin;//MAXV
int Max_Vtx,  Max_Iter; 
int f;
int fbest;

int *cruset;//MAXV
int len;
int tm1;
int tm2;
int *C0; //MAXV// the candidates node for ADD operation?
int *C1; //MAXV// the Candidate nodes for SWAP operation?
int *We; //MAXV// Weights of nodes?
int *BC;//MAXV // 
int len0; // the length of C0
int len1; // the length of C1
int *TC1;//MAXV // Temporal candidate nodes?
int Iter; // the number of iterations taken
int TABUL = 7;
int Wf;
int Wbest;
int *FC1;//MAXV
int *Tbest;//MAXV
int *TTbest;//MAXV
int Waim;
int Titer;

int len_best = 0;
int len_W;
int iter_best;
int Iteration[ 100 ];
double time_used[ 100 ];
int len_used[ 100 ];
int W_used[ 100 ];
char outfilename[30];
int len_improve;
int len_time;
int Wmode;
//int TABUL0 = 5;
//int *connect_clique_degree;
#ifdef clique_hash_mode
	CliqueHash *ptr_to_hashed_clique;
	int last_step_improved = 1;
#endif

#ifdef visited_vertex_set_mode
	Bijection *ptr_to_visited_vertex_set;
	//string last_op;
#endif

#ifdef aspiration_mode
	int aspired = 0;
#endif

#ifdef promising_starting_neighborhood_mode
	int biggest_found_clique_size = 0;
	int is_first_random_construction = 1;
#endif

#ifdef strategy_analysis_mode
	int restart_num = 0;
#endif

/************************************************************************/
/*   WYY: Variables for Configuration Checking                    */
/************************************************************************/
/* neighbor[i][j] = n means that i and n are connceted, and n is the jth 
 *  neighbor of n.
 */

/* time_stamp[m]=j means the conf_change value of node m recently 
 * changes at jth iteration.
 */

//struct  timeval start;
//struct  timeval end;

const long double EPS = 0.00001;

inline int long_double_equals(long double a, long double b, long double epsilon = EPS)
{
	return fabs(a - b) < epsilon;
}

inline int long_double_greater(long double a, long double b, long double epsilon = EPS)
{
	return a - b >= epsilon;
}

inline int long_double_smaller(long double a, long double b, long double epsilon = EPS)
{
	return b - a >= epsilon;
}

inline int long_double_geq(long double a, long double b, long double epsilon = EPS)
{
	return  (a - b >= epsilon || fabs(a - b) < epsilon);
}

inline int long_double_leq(long double a, long double b, long double epsilon = EPS)
{
	return (b - a >= epsilon || fabs(a - b) < epsilon);
}

void show_state()
{
	cout << "step: " << Iter << endl;	
	cout << "C0: " << endl;
	for(int i = 0; i < len0; i++) cout << C0[i] + 1 << '\t'; cout << endl;
	cout << "conf_change: " << endl;
	for(int i = 0; i < len0; i++) cout << conf_change[C0[i]] << '\t'; cout << endl;
	cout << "C1: " << endl;
	for(int i = 0; i < len1; i++) cout << C1[i] + 1 << '\t'; cout << endl;
	cout << "conf_change: " << endl;
	for(int i = 0; i < len1; i++) cout << conf_change[C1[i]] << '\t'; cout << endl;
	cout << "cand_solution: " << endl;
	for(int i = 0; i < len; i++) cout << cruset[i] + 1 << '\t'; cout << endl;

/*
	cout << "weights: " << endl;
	for(int i = 0; i < Max_Vtx; i++) cout << We[i] << '\t'; cout << endl;
*/
	cout << "clique_weight: " << Wf << endl;
#ifdef promising_starting_neighborhood_mode
	cout << "biggest_found_clique_size: " << biggest_found_clique_size << endl;
#endif
	cout << "========================================" << endl << endl;
	cout << "vertices: " << endl;
	for(int i = 0; i < Max_Vtx; i++) cout << i + 1 << '\t'; cout << endl;
	cout << "conf_change: " << endl;
	for(int i = 0; i < Max_Vtx; i++) cout << conf_change[i] << '\t'; cout << endl;
#ifdef clique_hash_mode
	cout << "last_step_improved: " << last_step_improved << endl;
#endif
#ifdef visited_vertex_set_mode
	cout << "visited_vertex_set: " << endl;
	for(int i = ptr_to_visited_vertex_set->begin(); i != ptr_to_visited_vertex_set->end(); i++)
		cout << ptr_to_visited_vertex_set->at(i) + 1 << '\t';
	cout << endl;
#endif
	cout << "*******************" << endl;
	getchar();
}

void dump_neighborhood();
#define DEBUG 0
///////////////////////////
inline int edge_is(int m, int n)
{
  int i;
  int index;
  for(i=0;i<neighbor_len[m];i++){
    index=neighbor[m][i];
    if(index==n)return 0;
  }
  return 1;
}

// section 0, initiaze
inline void Initializing()
{
     ifstream FIC;
     FILE *fp;
     FIC.open(File_Name);
     if ( FIC.fail() )
     {
           cout << "### Erreur open, File_Name " << File_Name << endl;
           getchar();
           exit(0);
     }
     char StrReading[100];
     //Max_Vclique=300;



     if ( FIC.eof() )
     {
           cout << "### Error open, File_Name " << File_Name << endl;
           exit(0);
     }
     int nb_vtx=0, nb_edg=-1, max_edg=0;
     int x1, x2;
	 int tt; 

    //FIC >> StrReading;


	do{
	FIC.getline(StrReading, 100);
	}while(StrReading[0] != 'p');
	char tmpstr1[10], tmpstr2[10];
	sscanf(StrReading, "%s %s %d %d", tmpstr1, tmpstr2, &Max_Vtx, &nb_edg);
//cout << Max_Vtx << " " << nb_edg << endl;
	//}while(tt != Max_Vtx);
	 //FIC >> tt >> Max_Vtx >> nb_edg;
	 //neighbor=(int **)malloc(nb_edg*2*sizeof(int));//neighbor set
	//neighbor = new int*[nb_edg * 2];
//cout << 1 << endl;
		neighbor = new int*[Max_Vtx];
//cout << 2 << endl;
/////////////////////////////////////
//cout << nb_edg << endl;
	edge = new Edge1[nb_edg];
//cout << 2.5 << endl;
	v_degree_tmp = new int[Max_Vtx];
	adjaclen = new int[Max_Vtx];//MAXV
//cout << 3 << endl;
	neighbor_len = new int[Max_Vtx];//MAXV
	conf_change = new int[Max_Vtx];//MAXV
	time_stamp = new int[Max_Vtx];//MAXV
//cout << 4 << endl;
	temp_array = new int[Max_Vtx];//MAXV
	vectex = new int[Max_Vtx];//MAXV
	funch = new int[Max_Vtx];//MAXV
//cout << 5 << endl;
	address = new int[Max_Vtx];//MAXV	
	tabuin = new int[Max_Vtx];//MAXV	
	cruset = new int[Max_Vtx];//MAXV	
	C0 = new int[Max_Vtx];//MAXV	
	C1 = new int[Max_Vtx];//MAXV	
	We = new int[Max_Vtx];//MAXV	
	BC = new int[Max_Vtx];//MAXV	
	TC1 = new int[Max_Vtx];//MAXV	
	FC1 = new int[Max_Vtx];//MAXV	
	Tbest = new int[Max_Vtx];//MAXV	
	TTbest = new int[Max_Vtx];//MAXV	
#ifdef clique_hash_mode
	ptr_to_hashed_clique = new CliqueHash(Max_Vtx);
#endif
	//connect_clique_degree = new int[Max_Vtx];
	//memset(connect_clique_degree, 0, sizeof(int)*Max_Vtx);
#ifdef visited_vertex_set_mode
	ptr_to_visited_vertex_set = new Bijection(Max_Vtx);
#endif
/////////////////////////////////////

	 for( int x = 0 ; x < Max_Vtx ; x++ ) 
	 {
		 conf_change[x] = 1; // initialize
         time_stamp[x] = 0;
		 adjaclen[x]=0;
         neighbor_len[x]=0;
         vectex[x]  = x;
         address[x] = x;
	}
	char tmpstr[10];
	int ii;
	for(ii = 0; ii < Max_Vtx; ii++)
	{
		int v_index;
		int v_weight;
		FIC >> tmpstr >> v_index >> v_weight;
		We[v_index - 1] = v_weight;
//cout << v_index << " " << v_weight << endl;
//getchar();
	}
	for(ii=0;ii<nb_edg;ii++)
	{


		FIC >> tmpstr >> x1 >> x2;

		x1--; x2--;
        if ( x1<0 || x2<0 || x1>=Max_Vtx || x2 >=Max_Vtx )
        {
            cout << "### Error of node : x1="
                 << x1 << ", x2=" << x2 << endl;
            exit(0);
        }
		max_edg++;
		neighbor_len[x1]++;
		neighbor_len[x2]++;
		edge[ii].v1 = x2;
		edge[ii].v2 = x1;     
		}
    int v;
    for (v=0; v<Max_Vtx; v++) adjaclen[v]=Max_Vtx-1-neighbor_len[v];//????

	for (v=0; v<Max_Vtx; v++) 
		//neighbor[v]=(int *)malloc( neighbor_len[v]*sizeof(int));
		neighbor[v]= new int[neighbor_len[v]];

    for(v=0; v<Max_Vtx; v++) v_degree_tmp[v]=0; 
	int e,v1,v2;  
	for (e=0; e<nb_edg; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;
		neighbor[v1][v_degree_tmp[v1]] = v2;
		neighbor[v2][v_degree_tmp[v2]] = v1;
		v_degree_tmp[v1]++;
		v_degree_tmp[v2]++;
	}
	int i;


     if ( 0 && max_edg != nb_edg )
     {
           cout << "### Error de lecture du graphe, nbre aretes : annonce="
                 << nb_edg << ", lu=" << max_edg  << endl;
           exit(0);
     }
     
     

     for( int x = 0; x < Max_Vtx; x++ )
     {
        //We[ x ] = (x+1)%Wmode + 1;
        BC[ x ] = 0;
        //We[ x ] = 1;
        //We[ x ] = ( rand() % 10 ) + 1;
     }
    
     FIC.close();


     
}

inline void free_memory()
{
	for(int v=0; v<Max_Vtx; v++)
		delete[] neighbor[v];
	delete[] neighbor;
	delete[] edge;
	delete[] v_degree_tmp;
	delete[] adjaclen;
	delete[] neighbor_len;
	delete[] conf_change;
	delete[] time_stamp;
	delete[]  temp_array;
	delete[] vectex;
	delete[]  funch;
	delete[] address;
	delete[] tabuin;
	delete[] cruset;
	delete[] C0;
	delete[] C1;
	delete[] We;
	delete[] BC;
	delete[] TC1;
	delete[] FC1;
	delete[] Tbest;
	delete[]  TTbest;
	//delete[] connect_clique_degree;
#ifdef clique_hash_mode
	delete ptr_to_hashed_clique;
#endif
#ifdef visited_vertex_set_mode
	delete ptr_to_visited_vertex_set;
#endif
}

// WYY
void dump_conf_change() {
	printf("\nconf_change:\n");
	for(int i = 0; i < Max_Vtx; i++) {
		printf("%4d(%d) ", i, conf_change[i]);
	}
	printf("\n");
}

// WYY
inline void neighbor_add(int node) {
	int node2;
	int num_neighbor = neighbor_len[node];
	
	//conf_change[node] = 0; // prevent this node from being removed immediately
	time_stamp[node] = Iter;
	for(int i = 0; i < num_neighbor; i++) {
		node2 = neighbor[node][i];
/*
		if(conf_change[node2] == 0){
			conf_change[node2] = 1;
			time_stamp[node2] = Iter;//????
		}
*/
		conf_change[node2] = 1;
	}
}
// WYY
inline void neighbor_drop(int node) {
	conf_change[node] = 0;
	time_stamp[node] = Iter;
}

// WYY
inline bool is_forbiden_cc(int node) {
	return (conf_change[node] == 0 ? true : false);
}


// WYY
void dump_neighborhood() {
	printf("Neighborhood:\n");
	for(int i = 0; i < Max_Vtx; i++) {
		printf(": ");
		for(int j = 0; j < neighbor_len[i]; j++)
			printf("%d ", neighbor[i][j]);
		printf("\n");	
	}
	return;
}

// WYY
void dump_cur_clique() {
	return;
	int n;
	printf("\ncurrent clique includes %d nodes:", len);
	for(int i = 0; i < len; i++) {
		n = cruset[i];
		printf("%d(%d) ", n, vectex[n]);
	}
	printf("\n");
}

inline int randomInt( int n )
{
    return rand() % n;
}

inline void backtract(int index);

inline void clearGamma()
{
    int i, j, k, l;
    memset( vectex, 0, tm1 );
    memset(  funch, 0, tm1 );
    memset(address, 0, tm1 );
    //memset( tabuin, 0, tm1 );
    for( i = 0; i < Max_Vtx; i++ )
    {
       C0[ i ] = i;// at the beginning all vertices can be added
       address[ i ] = i;
    }
    len0 = Max_Vtx;
    len1 = 0;
    len = 0;
    Wf = 0;
    //Wbest = 0;
#ifdef clique_hash_mode
		ptr_to_hashed_clique->reset_hash_entry();
#endif

}

inline void drop_clear()
{
	while(len != 0) backtract(0);
}

// WYY: C0 is the set of nodes that can be added? and C1 is the set nodes that can be swapped with?
int selectC0( ) // select a random vertex from add_set s.t. confChange=1, if called right after intialization, the confChange info is ignored
//breaking ties randomly
{
    //int i, j, k, l, m;
    //l = 0;
/*    
    if( len0 > 30 )
    {
       k = randomInt( len0 );
       return k;
    }
*/    
    // WYY: TC1 records the set of nodes which are not being forbidden
/*
    for( i = 0; i < len0; i++ )
    {
       k = C0[ i ];
       if( !is_forbiden_cc(k) ) // Added by: WYY
       {
         TC1[ l++ ] = i;
       }
    }
    
    if( l == 0 )
      return -1;
    else
    {
        k = randomInt( l );
        k = TC1[ k ];
        return k;
    }
*/
//cout << "in selectC0" << endl;
//for(int i = 0; i < len0; i++) cout << C0[i] << '\t'; cout << endl;
			if(len0 != 0)
				return rand() % len0;
			else
				return -1;
			
}
// WYY: Select from C0, a node, which is not in tabu and has the max weight, 
// or satisfy the aspiration rule though in tabu.
inline int WselectC0( )//breaking ties in favor of the greatest age
{
    int i, j, k, l1, l2, w1, w2, m;
    l1 = 0;
    l2 = 0;
    w1 = 0;
    w2 = 0;
    
    for( i = 0; i < len0; i++ )
    {
       k = C0[ i ];
       
       // WYY:store nodes that are not in tabu list and with the maximum weight, in FC1
       //if( !is_forbiden_cc(k) ) // Added by WYY
#ifdef aspiration_mode
				if(conf_change[k] || aspired == 1)
#else
				if(conf_change[k])
#endif
       {
           if( We[ k ] > w1 )
           {
              l1 = 0;
              w1 = We[ k ];
              FC1[ l1++ ] = i;
           }
           else if ( We[ k ] >= w1 )
           {
              FC1[ l1++ ] = i;// free best
           }
       }
/*
       else
       {  // WYY: stores nodes that are being in tabu but with the maximum weight, in TC1
           if( We[ k ] > w2 )
           {
              l2 = 0;
              w2 = We[ k ];
              TC1[ l2++ ] = i;
           }
           else if ( We[ k ] >= w2 )
           {
              TC1[ l2++ ] = i;// tabu best
           }
       }
*/
    }
#ifdef aspriation_mode
		aspired = 0;
#endif    
    // WYY: to check first if the aspiration rule is applicable.
    // If not, select a nodes which have the highest weithgts; break ties randomly.
#if 0
    if( (l2 > 0) && ( w2 > w1 ) && ((w2+Wf)>Wbest) ) //tabu best can lead to a new best clique
    {
        /*
        k = randomInt( l2 );
        k = TC1[ k ];
        */
        
        // WYY: Select the node with the oldest age
        k = TC1[0];
        int oldest_time = time_stamp[ C0[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l2; j++) {
			index = TC1[j];
			node = C0[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
			}
		}
		      
        //cout << "yes in aspiration w2+Wf = " << w2+Wf << endl;
        //getchar();
        return k;
    }  
    else 
#endif
		if( l1 > 0 ) //select free best
    {
    	/*
        k = randomInt( l1 );
        k = FC1[ k ];
        */
        
        // WYY
        k = FC1[0];
        int oldest_time = time_stamp[ C0[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l1; j++) {
			index = FC1[j];
			node = C0[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
				//cout << "elder one found" << endl;
				//getchar();
			}
		}
        return k;
    }
    else
    {
        return -1;
    }
}

// SelN: the index of the node selected from C0
inline void expand(int SelN)
{
    int i, j, k, k1, l, am, m, n, n1;
//cout << "in expanding func, selN: " << SelN << endl;    
    m = C0[ SelN ]; // the node is m
//cout << "in expanding func, to add " << m + 1 << endl; 
    cruset[ len++ ] = m; // add m into the current set
    vectex[ m ] = 1; // set the flag?
    Wf = Wf + We[ m ]; // Wf is the weight of the current clique, i,e, weight found. update it.

	 /* WYY: set the nodes in cruset as neighbors of m, and vice versa */
/*
	if(DEBUG) {
		printf("\nin expand");
		printf("\nadd node %d", m);
	}
*/
	neighbor_add(m);
	// remove from add-set    
    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;
    
    
    
    int node2;	
    for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    
    	//conf_change[m] = 0; // keep this node from being removed immediately
	time_stamp[m] = Iter;
    for(i=0;i<neighbor_len[m];i++){
     //node2 = neighbor[m][i];
/*
     if(conf_change[node2] == 0){
		conf_change[node2] = 1;
		time_stamp[node2] = Iter;// time stamp records at what time configration changes, different from the one described in the paper
	}
*/
			//conf_change[node2] = 1;
      n=neighbor[m][i];
      temp_array[n]=1;// temp_array[n]=1 means neighbors of the added vertex
    }

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m)continue;
       if(temp_array[i]==1)continue;
       n=i;
       funch[ n ]++; // WYY: funch[n] traces the number of nodes that are in the current clique
       
       if( funch[ n ] == 1 )
       {   // WYY: remove n from C0 (add-set)
           k1 = address[ n ];
           len0--;
           n1 = C0[ len0 ];
           C0[ k1 ] = n1;
           address[ n1 ] = k1;
           
           // put it into C1 (swap-set)
           C1[ len1 ] = n;
           address[ n ] = len1;
           len1++;
           
           BC[ n ] = m; // WYY: BC[n] = m denotes that n is Being Connected by m.// n is being paired with m
       }
       else if( funch[ n ] == 2 )
       {
           // remove n it from C1
           len1--;
           n1 = C1[ len1 ];
           k1 = address[ n ];
           C1[ k1 ] = n1;
           address[ n1 ] = k1;
       }
    }
#ifdef debug_mode
cout << "add " << m + 1 << endl;
#endif
//getchar();
#ifdef clique_hash_mode
		ptr_to_hashed_clique->update_hash_wrt_add(m);
#endif
#ifdef visited_vertex_set_mode
	if(!ptr_to_visited_vertex_set->element_in(m)) ptr_to_visited_vertex_set->insert_element(m);
#endif
//show_state();
#if 0    
    if( Wf > Wbest )
     {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
//cout << "Wbest: " << Wbest << endl;
		if(Wbest > lbest)
			iter_best = M_iter + Iter;
        /*
        for( i = 0; i < Max_Vtx; i++ )
        {
            Tbest[ i ] = vectex[ i ];
        }
        */
     }
#endif    
    //return 1;   
}

int selectC1( )// randomly select a confChanged vertex
{
    int i, j, k, l, m;
    l = 0;
    for( i = 0; i < len1; i++ )
    {
       k = C1[ i ];
       //if( !is_forbiden_cc(k) ) // WYY
				if(!conf_change[k])
       {
         TC1[ l++ ] = i;
       }
    }
    if( l == 0 )
      return -1;
    else
    {
        k = randomInt( l );
        k = TC1[ k ];
        return k;
    }
}

int kkk;
inline int WselectC1( )
{
     int i, j, k, l, l1, l2, wmn, w1, w2, m, n;
     l1 = 0;
     l2 = 0;
     w1 = -1000000;
     w2 = -1000000;
     l = 0;
     for( i = 0; i < len1; i++ )// for each pair in the swap-set
     {
         m = C1[ i ];
         n = BC[ m ];
         if( (vectex[ n ] == 1) && (edge_is( m, n)==1) )// vectex[n]==1 means in clique
           l++;
         else
         {
             for( j = 0; j < len; j++ )
             {
                k = cruset[ j ];
                if( edge_is(m, k)== 1 )
                  break;
             }
             BC[ m ] = k;
         }
     }//maintain BC
  //wyy-20150601
/*
     kkk=BMP;
	 if(len1<=BMP) kkk=len1;
*/	 
	 if(len1==0)return -1;
     //---add end
     //for( count = 0; count < kkk; count++ )
			for( i = 0; i < len1; i++ )
     {
		 
         m = C1[ i ];
         n = BC[ m ];
         wmn = We[ m ] - We[ n ];
		 //if( !is_forbiden_cc(m) ) // WYY//tabu best
#ifdef aspiration_mode
			//if(conf_change[m] || aspired == 1)
			if(conf_change[m])
#else
			if(conf_change[m])
#endif
         { // find the nodes that lead highest weight-increase for the current clique.
             if( wmn > w1 )
             {
                l1 = 0;
                w1 = wmn;
                FC1[ l1++ ] = i;
             }
             else if ( wmn >= w1 )
             {
                FC1[ l1++ ] = i;
             }
         }
/*
         else//free best
         {
             if( wmn > w2 )
             {
                l2 = 0;
                w2 = wmn;
                TC1[ l2++ ] = i;
             }
             else if ( wmn >= w2 )
             {
                TC1[ l2++ ] = i;
             }
         }
*/
     }
#ifdef aspiration_mode
		aspired = 0;
#endif
#if 0     
     if( (l2 > 0) && ( w2 > w1 ) && ((w2+Wf)>Wbest) )//free best
     {
        /* 
        k = randomInt( l2 );
        k = TC1[ k ];
        */
        
        // WYY: Select the oldest node
        k = TC1[0];
        int oldest_time = time_stamp[ C1[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l2; j++) {
			index = TC1[j];
			node = C1[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
			}
		}
        return k;
     }  
     else 
#endif
		if( l1 > 0 )//tabu best
     {
     /*
        k = randomInt( l1 );
        k = FC1[ k ];
       */  
        
        // WYY: Select the oldest node
        k = FC1[0];
        int oldest_time = time_stamp[ C1[k] ];
        int index;
        int node;
       	int time;
		for(int j = 1; j < l1; j++) {
			index = FC1[j];
			node = C1[index];
			time = time_stamp[ node ];
			if(time < oldest_time) {
				oldest_time = time;
				k = index;
				//cout << "elder nodes found " << endl;
				//getchar();
			}
		}
       	
       	return k;
     }
     else
     {
         return -1;
     }
}

inline int RselectC1()
{
		if(len1 > 0)
			return rand() % len1;
		else
			return -1;
}

inline void plateau( int SelN )
{
     int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;
     
     m = C1[ SelN  ];
     // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
     for(ti = 0; ti < len; ti++)
     {
         m1 = cruset[ ti ];
         if( edge_is(m1, m)== 1 )
            break;
     }
     
     Wf = Wf + We[ m ] - We[ m1 ];
     
     //the expand process, put m into the current independent set
     vectex[ m ] = 1;
     cruset[ len++ ] = m;

	 /* WYY: set the nodes in cruset as neighbors of m, and vice versa */
/*
	 if(DEBUG) {
	 	printf("\nin plateau: add node %d", m);
	 	dump_cur_clique();
	 }
*/
	 //neighbor_add(m); // Attention: here, we don't change conf_change values for m's neighbors
	 
     //delete m from C1
     k1 = address[ m ];
     len1--;
     n1 = C1[ len1 ];
     C1[ k1 ] = n1;
     address[ n1 ] = k1;
     
     
    for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    for(i=0;i<neighbor_len[m];i++){
     
      n=neighbor[m][i];
      temp_array[n]=1;
    }

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m)continue;
       if(temp_array[i]==1)continue;
       n=i;
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )//vectex[n]=0 means that n is not in clique
        {
             //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
             k1 = address[ n ];
             len0--;
             n1 = C0[ len0 ];
             C0[ k1 ] = n1;
             address[ n1 ] = k1;//remove from add-set
             
             C1[ len1 ] = n;
             address[ n ] = len1;
             len1++;
             BC[ n ] = m;//add into swap-set
           
             //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;//remove from swap-set
        }        
     } 
     
     //the backtrack process, delete m1 from the current independent set
     vectex[ m1 ] = 0;

	 len--;
     cruset[ ti ] = cruset[ len ];
     C1[ len1 ] = m1;
     address[ m1 ] = len1;
     len1++;
     
     /* WYY: neighborhood updating */
/*
	 if(DEBUG) {
	 	printf("\nin plateau: remove node %d", m1);
	 	dump_neighborhood();
	 }
*/
	 neighbor_drop(m1);
	 
    for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    for(i=0;i<neighbor_len[m1];i++){
     
      n=neighbor[m1][i];
      temp_array[n]=1;
    }//compute the neighbor set to exclude

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m1)continue;
       if(temp_array[i]==1)continue;
       n=i;
	 
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
           k1 = address[ n ];           
           len1--;
           n1 = C1[ len1 ];
           C1[ k1 ] = n1;
           address[ n1 ] = k1;//remove from swap-set
           
           C0[ len0 ] = n;
           address[ n ] = len0;
           len0++;//add into add-set
        }
        else if( funch[ n ] == 1 )//add into swap-set
        {
           C1[ len1 ] = n;
           address[ n ] = len1;
           len1++;
        }
     }
#ifdef debug_mode
cout << "add " << m + 1 << " and drop " << m1 + 1 << endl; 
#endif
//getchar();
#ifdef clique_hash_mode
		ptr_to_hashed_clique->update_hash_wrt_add(m);
		ptr_to_hashed_clique->update_hash_wrt_drop(m1);
#endif
#ifdef dynamic_tabu_management_mode
	if(!ptr_to_visited_vertex_set->element_in(m)) neighbor_add(m);
#endif
#ifdef aspiration_mode
	if(!ptr_to_visited_vertex_set->element_in(m)) aspired = 1;
#endif
#ifdef visited_vertex_set_mode
	if(!ptr_to_visited_vertex_set->element_in(m))
	{
#ifdef debug_mode
		cout << "new vertex: " << m + 1 << endl;
#endif 
		ptr_to_visited_vertex_set->insert_element(m);
	}
#endif
#if 0     
     if( Wf > Wbest )
     {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
//cout << "Wbest: " << Wbest << endl;
		if(Wbest > lbest)
			iter_best = M_iter + Iter;
        /*for( i = 0; i < Max_Vtx; i++ )
        {
            Tbest[ i ] = vectex[ i ];
        }*/
     }
#endif     
     
     //return 1;   
}

// WYY: find nodes with minimum weight in the current clique to remove
inline int Mumi_Weigt()
{
    int i, j, k, l1, m;
    int w1 = 5000000;
    l1 = 0;
    // WYY: find in cruset the nodes with lowest weights, breaking ties in favor of the oldest one
    for( i = 0; i < len; i++ )
    {
       k = cruset[ i ];
       if( We[ k ] < w1 )
       {
          l1 = 0;
          w1 = We[ k ];
          FC1[ l1++ ] = i;
       }
       else if ( We[ k ] <= w1 )
       {
          FC1[ l1++ ] = i;
       }
    }

#ifdef aspiration_mode
		aspired = 0;
#endif
    
    if( l1 == 0 )
      return -1;
    /*
    k = randomInt( l1 );
    k = FC1[ k ];
    */
    
    // WYY: remove an oldest node
    k = FC1[0];
    int oldest_time = time_stamp[ cruset[k] ];
    int cur_index;
    int cur_node;
    int cur_time;
    for(int i = 1; i < l1; i++) {
    	cur_index = FC1[i];
    	cur_node = cruset[cur_index ];
    	cur_time = time_stamp[cur_node];
    	if(cur_time < oldest_time) {
    		oldest_time = cur_time;
    		k = cur_index;
    		//cout << "elder node found " << endl;
    		//getchar();
    	}
    }
    return k;
}

inline int rand_index_in_clique()
{
		if(len) return rand() % len;
		return -1;
}

inline void backtract(int SelN)
{
     int i, j, k, l, m, m1, n, k1, n1;
     //int ti = Mumi_Weigt();
     //if( ti == -1 )
      //return -1;
      
     //if(DEBUG) printf("in backtrack");
     m1 = cruset[SelN];
     //m1 = cruset[ ti ];
     Wf = Wf - We[ m1 ];
     vectex[ m1 ] = 0;

     len--;
     //cruset[ ti ] = cruset[ len ];
			cruset[SelN] = cruset[len];
     
     /* WYY: functions of neighborhood updating */
	 neighbor_drop(m1);
/*
	 if(DEBUG) {
	 	printf("\nremove node %d", m1);
	 	dump_cur_clique();
	 }
*/
	 //dump_conf_change();
	 //getchar();
     
     C0[ len0 ] = m1;
     address[ m1 ] = len0;
     len0++;
     
   for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
    for(i=0;i<neighbor_len[m1];i++){
     
      n=neighbor[m1][i];
      temp_array[n]=1;
    }

    for(i=0;i<Max_Vtx;i++)
    {
       if(i==m1)continue;
       if(temp_array[i]==1)continue;
       n=i;     

        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
           k1 = address[ n ];           
           len1--;
           n1 = C1[ len1 ];
           C1[ k1 ] = n1;
           address[ n1 ] = k1;//remove from swap-set
           
           C0[ len0 ] = n;
           address[ n ] = len0;
           len0++;//add into add-set
        }
        else if( funch[ n ] == 1 )
        {
           C1[ len1 ] = n;
           address[ n ] = len1;
           len1++;// add into add-set
        }
     }
#ifdef debug_mode
cout << "drop " << m1 + 1 << endl;
#endif
//getchar();
#ifdef clique_hash_mode
		ptr_to_hashed_clique->update_hash_wrt_drop(m1);
#endif
			//return 1;
}

#ifdef promising_starting_neighborhood_mode
inline int clique_weight_estimation(int v)
{
	int sum_weight = We[v];
	for(int i = 1; i < biggest_found_clique_size; i++)//this size wierd, too restricted
	{
		sum_weight += We[neighbor[v][rand() % neighbor_len[v]]];
	}
cout << "considering vertex: " << v + 1 << endl;
cout << "estimated_weight: " << sum_weight << endl;
	return sum_weight;
}
#endif

inline void random_construction(int &l)
{
#ifdef visited_vertex_set_mode
		ptr_to_visited_vertex_set->clear();
		//last_op = "add";
#endif
		int am;

#ifdef promising_starting_neighborhood_mode
		if(!is_first_random_construction)
		{
			while(1)
			{
				am = selectC0();
				int estimated_weight = clique_weight_estimation(C0[am]);
				if(double(rand() % 10000) / 10000.0 < double(estimated_weight) / double(lbest))
				{
					expand(am);
					Iter++;
					break;
				}
			}
		}
		is_first_random_construction = 0;
#endif
			
     while( 1 )
     {
        am = selectC0();
        if( am != -1 )
        {
            expand( am );
            Iter++;
            //if( Wbest == Waim )
               //return Wbest;
        }
        else 
            break;
     }
#ifdef clique_hash_mode
			last_step_improved = 1;
#endif    
//show_state();
//getchar();
		if(Wf > Wbest)
		{
			Wbest = Wf;
			len_W = len;
			for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
			iter_best = Iter;
     times(&finish);
     double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
	 finish_time = round(finish_time * 100)/100.0;
			real_solve1 = finish_time;
		}
#if 0
	 if(Wbest>lbest){lbest=Wbest; len_W = len; real_solve1 = real_solve2; iter_best = M_iter + Iter;}
	 if(finish_time>time_limit){M_iter+=Iter;

	 //printf("%.2f	%d %d\n", real_solve1,lbest, M_iter);
			cout << "o " << lbest << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(M_iter) / 1000.0 / finish_time << endl;	
	 exit(0);}
#endif
}

inline int tabu( int Max_Iter )
{
     int i, j, k, l, bestlen = 0, am, am1, ww, ww1, ww2, ti, m1;
     Iter = 1;
     clearGamma();
//cout << "initializing" << endl;
//show_state();
#ifdef clique_hash_mode
		last_step_improved = 1; 
#endif
/*
     while( 1 )
     {
        am = selectC0();
        if( am != -1 )
        {
            l = expand( am );
            Iter++;
            if( Wbest == Waim )
               return Wbest;
        }
        else 
            break;
     }
     times(&finish);
     double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
	 finish_time = round(finish_time * 100)/100.0;
	 if(Wbest>lbest){lbest=Wbest; len_W = len; real_solve1 = real_solve2; iter_best = M_iter + Iter;}
	 if(finish_time>time_limit){M_iter+=Iter;

	 //printf("%.2f	%d %d\n", real_solve1,lbest, M_iter);
			cout << "o " << lbest << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(M_iter) / 1000.0 / finish_time << endl;	
	 exit(0);}
*/     
			//random_construction(l);
  
     //while( Iter < Max_Iter )
			while(1)
     {
				if(Iter % 1000 == 0)
				{	
     			times(&finish);
     			double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
	 				finish_time = round(finish_time * 100)/100.0;
					if(finish_time > time_limit) break;					
				}
				if(len == 0)
				{
					random_construction(l);
#ifdef clique_hash_mode
					last_step_improved = 1;
#endif			
				}
        am = WselectC0();
				
				//if(0)//anyway, compare all possible ways of adding and swapping, then choose the score-greatest way
				//if(len == 1 && am != -1)//when |C|=1 and add-set is not empty, forbid swapping, base-2
				if(len == 1)//when |C|=1, forbid swapping, base-1
				{
					am1 = -1;
				}
				else
				{
        	am1 = WselectC1();
				}

#ifdef first_improved_revisit_restart_mode
				int is_first_improved_step = 0;
#endif

		    ww = We[ C0[ am ] ];
		    ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
//should test whether restart is needed
		      
				if(am != -1)
				{
#ifdef first_improved_revisit_restart_mode
					if(!last_step_improved)
					{ 
						is_first_improved_step = 1;
#ifdef debug_mode
						cout << "is_first_improved_step" << endl;
#endif
					}
#endif
		      if( am1 != -1 )
		      {
		          if( ww > ww1 )
		          {
		              expand( am );  
		              Iter++;
		          }
		          else
		          {
		              plateau( am1 );
		              Iter++;
		          }
		      }
		      else if( am1 == -1 )//add-set is not empty; swap-set is empty
		      {
		           expand( am );
		              
		           Iter++;
		      }
#ifdef clique_hash_mode
					last_step_improved = 1;
#endif
				}
        else
				{
#ifdef promising_starting_neighborhood_mode
					if(biggest_found_clique_size < len) biggest_found_clique_size = len;
#endif
#ifdef clique_hash_mode
					if(am1 == -1 || ww1 < 0)
					{
	#ifdef local_optimum_revisit_restart_mode
						if(last_step_improved)
						{
							ptr_to_hashed_clique->mark_hash_entry();	
	//cout << "mark hash" << endl;
							if(ptr_to_hashed_clique->curr_hash_entry_marked_too_frequently())
							{
	//cout << "hash bound exceeds" << endl;
#ifdef strategy_analysis_mode
restart_num++;
#endif
								//clearGamma();
								drop_clear();
								Iter++;
								continue;

							}
						}
	#endif
						last_step_improved = 0;
					}
					else
					{
	#ifdef first_improved_revisit_restart_mode
						if(!last_step_improved)
						{
							is_first_improved_step = 1;
#ifdef debug_mode
							cout << "is_first_improved_step" << endl;
#endif
						}
	#endif
						last_step_improved = 1;
					}					
#endif	
		      ti = Mumi_Weigt();					 
					//if( (am == -1) && (am1 != -1) )//add-set is empty; swap-set is not empty
					if( am1 != -1 )//add-set is empty; swap-set is not empty
		      {
		           m1 = cruset[ ti ];
		           //ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
		           ww2 = - We[ m1 ];
		           if( ww1 > ww2 )
		           {
		              plateau( am1 );
		              Iter++;
		           }
		           else
		           {
		               backtract(ti);		               
		               Iter++;
		           }
		      }
		      else//add-set is empty; swap-set is empty
		      {
#ifdef rand_drop_mode
							 backtract(rand_index_in_clique());
#else
		           backtract(ti); 
#endif   
		           Iter++;
		      }
				}
#ifdef debug_mode
show_state();
#endif

#ifdef first_improved_revisit_restart_mode
				if(is_first_improved_step)
				{
					ptr_to_hashed_clique->mark_hash_entry();
					if(ptr_to_hashed_clique->curr_hash_entry_marked_too_frequently())
					{
#ifdef debug_mode
cout << "this solution is revisited" << endl;
#endif
#ifdef strategy_analysis_mode
restart_num++;
#endif
						drop_clear();
						Iter++;
						continue;
					}
				}
#endif
				if(Wf > Wbest)
				{
					Wbest = Wf;
					len_W = len;
					for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
					iter_best = Iter;
        times(&finish);
        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
		finish_time = round(finish_time * 100)/100.0;
				real_solve1 = finish_time;
				}
#if 0
        //times(&finish);
        //double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
		//finish_time = round(finish_time * 100)/100.0;
        if(finish_time > time_limit){
			M_iter+=Iter;
        	if(Wbest>lbest){
					lbest=Wbest;
					real_solve1 = real_solve2;
				}
			//printf("%.2f	%d %d\n", real_solve1, lbest, M_iter);
			cout << "o " << lbest << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(M_iter) / 1000.0 / finish_time << endl;
			
			exit(0);
		}
#endif
     }

     return Wbest;
}

void verify()
{
     int i, j, k1, k2, l, m;
     for( i = 0; i < Max_Vtx; i++ )
     {
          if( TTbest[ i ] == 1 )
          {
              for( j = i+1; j < Max_Vtx; j++ )
              if( (TTbest[ j ] == 1) && ( edge_is(i, j)== 1 ) )
                  cout << "hello there is something wrong" << endl;
          }
     }
     cout << "verified" << endl;
}

// WYY: Validate that the cruset is indeed a clique
void validate() {
	int i, j, k1, k2, l, m;
     for( i = 0; i < Max_Vtx; i++ )
     {
          if( vectex[ i ] == 1 )
          {
              for( j = i + 1; j < Max_Vtx; j++ )
              if( (vectex[ j ] == 1) && ( edge_is(i, j)== 1 ) ) {
                  cout << "hello there is something wrong" << endl;
              	getchar();
              }
          }
     }
}

void Output()
{
    int i , j, k, l, sum; 
    FILE *fp ;
    int len = strlen(File_Name);
    
    strcpy(outfilename,File_Name) ;
    outfilename[len]='.';
    outfilename[len+1]='o';
    outfilename[len+2]='u';
    outfilename[len+3]='t';
    outfilename[len+4]='\0';

    fp = fopen(outfilename, "a+"); 
    for( i = 0; i < 1; i++ )
    {
        fprintf(fp, "sum = %6d, iter = %6d, len = %5d,  time = %8lf \n", 
        	W_used[ i ], Iteration[ i ], len_used[ i ],  time_used[ i ] ); 
    }
    
    fclose(fp); // WYY
    return;
    
    fprintf(fp, "\n\n the total information: \n");
    int wavg, iteravg, lenbb, success;
    wavg = iteravg = lenbb = success = 0;
    int best_v = 0;
    double timeavg = 0; 
    for( i = 0; i < 100; i++ )
		if( W_used[ i ] > best_v )
		{
			best_v = W_used[ i ];  
			lenbb = len_used[ i ];
		}
    
    int count = 0;
    fprintf(fp, "\n The best weight value for the maximum weighted problem is %6d \n", best_v);
    for( i = 0; i < 100; i++ )
    {
       wavg = wavg + W_used[ i ];
    }  
    double twavg = (double (wavg)) / 100 ; 
    for( i = 0; i < 100; i++ )
    if( W_used[ i ] == best_v )
    {
        count++;
        iteravg = iteravg + Iteration[ i ];
        timeavg = timeavg + time_used[ i ];
    }
    
    iteravg =  int ( (double (iteravg)) / count );
    timeavg = timeavg / (count*1000);
    fprintf(fp, "avg_sum = %10lf, succes = %6d, len = %5d, avg_iter = %6d,  time = %8lf \n", 
    			twavg, count, lenbb,  iteravg, timeavg );
    fclose(fp) ;
    return ;
}

inline void Max_Tabu()
{
     int i, j, k, l, m;
     //lbest = 0;
     //int lenbest = 0;
     //Titer = 0;
	 times(&start);
////////////////////////////////

/////////////////////////////////
         l = tabu(len_improve);
/*
         M_iter = M_iter + Iter; 
         if( l > lbest )
         {
			 real_solve1=real_solve2;
			 lbest = l; 
			 Titer = M_iter; 
			 len_W = len_best;        
         }
         
         if( l >= Waim )
           return lbest;
*/
         //cout << " l = " << l << " i = " << i << endl;
		 
		 // wyy: clear configuration Information for the next restart
/*
		 for(int j = 0; j < Max_Vtx; j++) {
          	conf_change[j] = 1;
          	time_stamp[j] = 0;
         }
*/
		 times(&finish);
		 double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
		 finish_time = round(finish_time * 100)/100.0;
		 if(finish_time>time_limit){
				//printf("%.2f	%d\n", real_solve1,lbest);
			cout << "o " << Wbest << endl;
			cout << "c size " << len_W << endl;
			cout << "c solveTime " << real_solve1 << endl;
			cout << "c searchSteps " << iter_best << endl;
			cout << "c stepSpeed(/ms) " << double(Iter) / 1000.0 / finish_time << endl;
#ifdef strategy_analysis_mode
			cout << "restart_num: " << restart_num << endl;
			cout << "restart freq: " << double(Iter) / restart_num << endl;
#endif
			}
 
}

int main(int argc, char **argv)
{
     int seed;
	 File_Name = argv[1];
	 //Wmode=200;
	 time_limit=atof(argv[3]);
	 //len_improve=4000;
	 seed = atoi(argv[2]);
	 //printf("%s	",argv[1]);
	 //Waim=INT_MAX;
	 srand(seed);
	 Initializing();
	 tm1 = Max_Vtx*sizeof( int );
     //  cout << "finish reading data" << endl;
     int i, l;
	 //len_time = (int (100000000 / len_improve) ) + 1;
    //   cout << "len_time = " << len_time << endl;
     Max_Tabu();
		verify();
		free_memory();
}
