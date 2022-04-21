#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
NumericVector SP(CharacterMatrix ss){
  int r = ss.nrow();
  int c = ss.ncol();
  NumericVector ans(r);
  if(c == 1){
    return ans;
  }
  for(int i=0;i<r;i++){
    ans(i) = 0;
    for(int k=0;k<c;k++){
      for(int j=0;j<c;j++){
        if(k>=j){
          continue;
        }
        ans(i) += ss(i,k)==ss(i,j) ? (ss(i,k)=='-'?0:5) : (ss(i,k)=='-'||ss(i,j)=='-'?-5:-4);
      }
    }
  }
  return ans;
}

// [[Rcpp::export]]
CharacterMatrix liang_NW(CharacterMatrix s1, CharacterMatrix s2,Function print){
  int lenr = s1.nrow() + 1; //单序列
  int lenc = s2.nrow() + 1; //单/多序列
  int numc = s2.ncol();
  NumericMatrix marks(lenr, lenc); //构建得分矩阵
  NumericMatrix ways(lenr, lenc); //回溯矩阵
  NumericVector sp;
  sp = SP(s2);
  int gap = -5;
  int now_mark,left,up,go;
  for(int i = 0;i<lenr;i++){
    for(int k = 0;k<lenc;k++){
      if(i==0 || k==0){
        if(i==0&&k==0){
          marks(i,k) = 0;
        } else if(i==0){
          marks(i,k) = marks(i,k-1) + gap*numc + sp(k-1);
          ways(i,k) = 3;
        } else {
          marks(i,k) = marks(i-1,k) + gap;
          ways(i,k) = 2;
        }
        continue;
      }
      now_mark = sp(k-1);
      for(int q=0;q<numc;q++){
        now_mark += s2(k-1,q)=="-" ? gap : (s1(i-1,0) == s2(k-1,q) ? 5 : -4);
      }
      left = marks(i,k-1) + gap;
      up = marks(i-1,k) + gap;
      go = marks(i-1,k-1) + now_mark;
      marks(i,k) = (left>up?left:up) > go ? (up>left ? up : left) : go;
      ways(i,k) = (left>up?left:up) > go ? (up>left ? 3 : 1) : 2;
    }
  }
  //回溯
  int nowi = lenr - 1;
  int nowk = lenc - 1;
  int way[lenr+lenc]; //1-left, 2-go, 3,up
  int far = 0;
  int addSpace = -1;
  while(true){
    if(nowi==0 || nowk==0){
      if(nowi==0 && nowk==0){
        break;
      }
      way[far] = nowi==0 ? 1 : 3;
    } else {
      //up = marks(nowi-1,nowk);
      //left = marks(nowi,nowk-1);
      //go = marks(nowi-1,nowk-1);
      //way[far] = (left>up?left:up) > go ? (up>left?3:1) : 2;
      way[far] = ways(nowi,nowk);
    }
    nowi = way[far]==1 ? nowi : nowi - 1;
    nowk = way[far]==3 ? nowk : nowk - 1;
    addSpace = way[far]==1 ? addSpace + 1 : addSpace;
    far += 1;
  }
  CharacterMatrix ans(far,numc+1);
  //print(far - addSpace - lenr);
  //return(ans);
  nowi = lenr-2;
  nowk = lenc-2;
  int now_far = far-1;
  for(int i=0;i<far;i++){
    if(way[i]==1){
      ans(now_far,0) = "-";
      
      for(int j=1;j<numc+1;j++){
        ans(now_far,j) = s2(nowk,j-1);
      }
      //nowk = nowk>0 ? nowk - 1:nowk;
      nowk -= 1;
    }
    if(way[i]==2){
      ans(now_far,0) = s1(nowi, 0);
      //nowi = nowi>0 ? nowi - 1:nowi;
      nowi -= 1;
      
      for(int j=1;j<numc+1;j++){
        ans(now_far,j) = s2(nowk,j-1);
      }
      //nowk = nowk>0 ? nowk - 1:nowk;
      nowk -= 1;
    }
    if(way[i]==3){
      ans(now_far,0) = s1(nowi, 0);
      //nowi = nowi>0 ? nowi - 1:nowi;
      nowi -= 1;
      
      for(int j=1;j<numc+1;j++){
        ans(now_far,j) = "-";
      }
    }
    now_far -= 1;
  }
  return ans;
}

// [[Rcpp::export]]
CharacterMatrix tranGen_c(CharacterMatrix tm, NumericMatrix gen){
  int nr = tm.nrow();
  int nc = tm.ncol();
  int ch,p,r,nowk,nowj;
  CharacterVector v(nr),tmp(nr);
  for(int i=0;i<nc;i++){
    ch = i;
    p = gen(0,i) - 1;
    r = gen(1,i);
    v = tm(_,ch);
    nowk = 0;
    nowj = 0;
    if(r==p) continue;
    for(int k=0;k<nr;k++)
    {
      if(k==p)
      {
        nowk++;
        tmp(nowj) = v(nowk);
        nowj++;
        nowk++;
        continue;
      }
      if(k!=r){
        tmp(nowj) = v(nowk);
        nowk++;
        nowj++;
        continue;
      }
      else
      {
        tmp(nowj) = v(p);
        nowj++; 
      }
    }
    tm(_,ch) = tmp;
  }
  return(tm);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
