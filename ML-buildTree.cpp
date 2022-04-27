#include <Rcpp.h>
#include <math.h>
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
double get_random(){
  return (double)(rand() % 100) / 100;
}

class Ind
{
public:
  Ind *father;
  Ind *childs[2];
  int leaf_id = -1;
  int child_num;
  int update_num(int id);
  int branch_long = 0;
  int base_long = 0;
  CharacterVector branch_base(int np);
  CharacterMatrix seqs;
  int front_ergodic(int id, Ind *hx[]);
  Ind *recommend();
  Ind *add_father(int leaf);
  void right_lerp(); //右倾
  void left_mesh(Ind *top, Function print); //向左合并
  void delete_father();
  Ind(int id, Ind *fa, CharacterMatrix ss);
  ~Ind();
};

int cal_tree(Ind *start){
  int nps = start->seqs.nrow();
  int tree_long = 0;
  for(int i=0;i<nps;i++){
    start->branch_base(i);
    tree_long += start->branch_long;
  }
  return tree_long;
}

int cal_tree2(Ind *start){
  int nps = start->seqs.nrow();
  int tree_long = 0;
  for(int i=0;i<nps;i++){
    start->branch_base(i);
    tree_long += start->base_long;
  }
  return tree_long;
}

Ind::Ind(int id, Ind *fa, CharacterMatrix ss){
  leaf_id = id;
  father = fa;
  seqs = ss;
}
Ind::~Ind(){
  if(leaf_id==-1){
    for(int i=0;i<2;i++){
      if(childs[i]!=NULL){
        delete childs[i]; 
      }
    } 
  }
}

void Ind::right_lerp(){
  if(leaf_id!=-1) return;
  childs[0]->right_lerp();
  childs[1]->right_lerp();
  if(childs[0]->child_num<=childs[1]->child_num) return;
  Ind *tmp = childs[0];
  childs[0] = childs[1];
  childs[1] = tmp;
}

void Ind::left_mesh(Ind *top, Function print){
  if(leaf_id != -1 || childs[1]->leaf_id!=-1) return;
  childs[0]->left_mesh(top,print);
  int mark = cal_tree(top);
  int mark2 = cal_tree2(top);
  Ind *fa = new Ind(-1, this, seqs);
  childs[0]->father = fa;
  fa->childs[0] = childs[0];
  childs[0] = fa;
  
  Ind *c_node = childs[1];
  childs[1]->childs[0]->father = fa;
  fa->childs[1] = childs[1]->childs[0];
  childs[1]->childs[0] = NULL;
  childs[1] = childs[1]->childs[1];
  c_node->childs[1] = NULL;
  delete c_node;
  
  top->update_num(0);
  int after_mark = cal_tree(top);
  if(after_mark<mark){
    fa->right_lerp();
    this->left_mesh(top,print);
  } else if(after_mark == mark && cal_tree2(top) - mark2 < 10) {
    //print(cal_tree2(top) - mark2);
    fa->right_lerp();
    this->left_mesh(top,print);
  } else {
    childs[0] = fa->childs[0];
    fa->childs[0] = NULL;
    childs[0]->father = this;
    
    Ind *a_node = new Ind(-1, this, seqs);
    a_node->childs[0] = fa->childs[1];
    fa->childs[1] = NULL;
    delete fa;
    a_node->childs[1] = childs[1];
    childs[1] = a_node;
    childs[1]->left_mesh(top,print);
  }
}

Ind* Ind::add_father(int leaf){
  Ind *tmp = father;
  father = new Ind(-1, tmp, seqs);
  if(tmp!=NULL){
    for(int i=0;i<2;i++){
      if(tmp->childs[i]==this){
        tmp->childs[i] = father;
        break;
      }
    }
  }
  father->childs[0] = this;
  father->childs[1] = new Ind(leaf, father, seqs);
  return(father);
}

void Ind::delete_father(){
  Ind *tmp = father->father;
  if(tmp!=NULL){
    for(int i=0;i<2;i++){
      if(tmp->childs[i]==father){
        tmp->childs[i] = this;
        break;
      }
    }
  }
  for(int i=0;i<2;i++){
    if(father->childs[i]==this){
      father->childs[i] = NULL;
      break;
    }
  }
  delete father;
  father = tmp;
}

int Ind::front_ergodic(int id, Ind *hx[]){
  if(leaf_id==-1){
    id = childs[0]->front_ergodic(id, hx);
    id = childs[1]->front_ergodic(id, hx);
  }
  hx[id] = this; 
  id++;
  return id;
}

int Ind::update_num(int num){
  int n1 = 0,n2 = 0;
  if(leaf_id==-1){
    n1 = childs[0]->update_num(num);
    n2 = childs[1]->update_num(num);
    child_num = n1+n2;
  } else {
    child_num = 0; 
  }
  num += n1 + n2 + 1;
  return num;
}

Ind* Ind::recommend(){
  double p = get_random();
  if(leaf_id==-1 && p < 1/(double)child_num){
    bool is_left = childs[0]->leaf_id==-1;
    bool is_right = childs[1]->leaf_id==-1;
    if(is_left && is_right){
      if(p > 0.5){
        return childs[0]->recommend();
      } else {
        return childs[1]->recommend();
      }
    } else {
      if(is_left){
        //右边是叶子
        return childs[0]->recommend();
      } else if(is_right){
        return childs[1]->recommend();
      } else {
        return this;
      }
    }
  }
  return this;
}


CharacterVector Ind::branch_base(int np){
  if(leaf_id!=-1){
    CharacterVector new_base(1);
    new_base[0] = seqs(np,leaf_id);
    branch_long = 0;
    base_long = 1;
    return new_base;
  }
  CharacterVector str1 = childs[0]->branch_base(np);
  CharacterVector str2 = childs[1]->branch_base(np);
  int len1 = str1.size();
  int len2 = str2.size();
  int sam_num = 0;
  for(int i=0;i<len1;i++){
    for(int k=0;k<len2;k++){
      if(str1[i]==str2[k]){
        sam_num++;
      }
    }
  }
  if(sam_num){
    branch_long = childs[0]->branch_long + childs[1]->branch_long;
    base_long = childs[0]->base_long + childs[1]->base_long + sam_num;
    CharacterVector new_base(sam_num);
    int nowi = 0;
    bool not_have = true;
    for(int i=0;i<len1;i++){
      for(int k=0;k<len2;k++){
        if(str1[i]==str2[k]){
          not_have = true;
          for(int j=0;j<nowi;j++){
            if(str1[i] == new_base[j]){
              not_have = false;
              break;
            }
          }
          if(not_have){
            new_base[nowi] = str1[i];
            nowi++; 
          }
        }
      }
    }
    return new_base;
  } else {
    branch_long = 1 + childs[0]->branch_long + childs[1]->branch_long;
    base_long = childs[0]->base_long + childs[1]->base_long + len1 + len2;
    CharacterVector new_base(len1+len2);
    int nowi = 0;
    for(int i=0;i<len1;i++){
      new_base[nowi] = str1[i];
      nowi++; 
    }
    for(int i=0;i<len2;i++){
      new_base[nowi] = str2[i];
      nowi++;
    }
    return new_base;
  }
}

bool change_node(Ind *child1, Ind *child2){
  Ind *fa1 = child1->father;
  Ind *fa2 = child2->father;
  if(fa1 == fa2){
    return false;
  }
  for(int i=0;i<2;i++){
    if(fa1->childs[i]==child1){
      fa1->childs[i] = child2;
      break;
    }
  }
  for(int i=0;i<2;i++){
    if(fa2->childs[i]==child2){
      fa2->childs[i] = child1;
      break;
    }
  }
  child1->father = fa2;
  child2->father = fa1;
  return true;
}

// [[Rcpp::export]]
NumericVector MP_build(CharacterMatrix seqs, Function print) {
  int nc = seqs.ncol();
  Ind *start = new Ind(-1, NULL, seqs);
  Ind *nowExtend = start;
  for(int i=0;i<nc-2;i++){
    nowExtend->childs[0] = new Ind(-1, nowExtend, seqs);
    nowExtend->childs[1] = new Ind(i, nowExtend, seqs);
    nowExtend = nowExtend->childs[0];
  }
  nowExtend->childs[0] = new Ind(nc-2, nowExtend, seqs);
  nowExtend->childs[1] = new Ind(nc-1, nowExtend, seqs);
  
  int apoint = nc + nc - 1;
  int mark,deta,num1,num2;
  Ind *back_list[apoint];
  Ind *nowP,*ind1,*ind2;
  
  double Tp = 100;
  double Te = 1;
  double alpha = 0.99;
  int cnt = 10;
  int bestMark = cal_tree(start);
  //print(bestMark);
  //double gt = 0.3;
  
  start->update_num(0);
  while(Tp > Te){
    cnt = 10;
    while(cnt > 0){
      cnt--;
      
      /***
      if(get_random()>0.5){
        start->front_ergodic(0,back_list);
        int nowi = 0;
        int one_not = 0;
        while(one_not < 2){
          nowi = 1;
          one_not = 0;
          for(int k=0;k<apoint;k++){
            if(back_list[k]->leaf_id!=-1){
              nowi++;
              if(get_random() < (double)nowi/(double)nc){
                if(one_not==0){
                  ind1 = back_list[k];
                  one_not = 1;
                } else {
                  ind2 = back_list[k];
                  one_not = 2;
                  break;
                }
              }
            }
          }
        }
      } else {
        nowP = start->recommend(gt);
        ind1 = nowP->childs[0]->recommend(0.01);
        ind2 = nowP->childs[1]->recommend(0.01);
      }
       ***/
      
      nowP = get_random() > 0.5 ? start->recommend() : start;
      if(get_random() > 0.5)
      {
        ind1 = nowP->childs[0]->recommend();
        ind2 = nowP->childs[1]->recommend();
      } else {
        num1 = nowP->childs[0]->child_num+1;
        num2 = nowP->childs[1]->child_num+1;
        Ind *l1[num1];
        Ind *l2[num2];
        nowP->childs[0]->front_ergodic(0, l1);
        nowP->childs[1]->front_ergodic(0, l2);
        ind1 = l1[rand() % num1];
        ind2 = l2[rand() % num2]; 
      }
      
      if(!change_node(ind1, ind2)){
        cnt++;
        continue;
      }
      mark = cal_tree(start);
      deta = bestMark - mark;
      // || exp((double)deta/Tp)
      if(deta > 0 || (get_random() < exp((double)deta/Tp) && Tp > 10)){
        bestMark = deta > 0 ? mark : bestMark;
        cnt = deta > 0 ? 10 : cnt;
        //bestMark = mark;
        //cnt = 10;
        start->update_num(0);
      } else {
        change_node(ind1, ind2);
      }
    }
    
    Tp = Tp * alpha;
  }
  
  print(bestMark);
  start->front_ergodic(0,back_list);
  NumericVector ans(apoint);
  for(int i =0;i<apoint;i++){
    ans[i] = back_list[i]->leaf_id;
  }
  delete start;
  return ans;
}

// [[Rcpp::export]]
NumericVector qian_zhui_he(NumericVector tree){
  int len = tree.size();
  NumericVector ans(len);
  int qian = 0;
  for(int i=0;i<len;i++){
    qian = tree[i]==0 ? qian + 1 : 0;
    ans[i] = qian;
  }
  return ans;
}

bool place_the_seq(Ind* grandfa,Ind* now_test, int id, int now_long,Function print){
  if(now_test->leaf_id != -1){
    return false;
  }
  Ind *left = now_test->childs[0];
  Ind *right = now_test->childs[1];
  left->add_father(id);
  int left_mark = cal_tree(now_test);
  //int left_mark = cal_tree(grandfa);
  left->delete_father();
  right->add_father(id);
  int right_mark = cal_tree(now_test);
  //int right_mark = cal_tree(grandfa);
  bool go_left = right_mark > left_mark ? true : false;
  if(go_left){
    left->add_father(id);
    right->delete_father();
  }
  int min_mark = cal_tree(grandfa);
  if(go_left){
    left->delete_father();
  } else {
    right->delete_father();
  }
  if(min_mark > now_long){
    return false;
  }
  bool not_best = true;
  if(go_left){
    not_best = place_the_seq(grandfa, left, id, min_mark, print);
  } else {
    not_best = place_the_seq(grandfa, right, id, min_mark, print);
  }
  if(!not_best){
    /*
    if(grandfa == now_test){
      //print('d');
      //print(id);
      //print(left_mark);
      //print(right_mark);
      //print('d');
      //return false;
    }
    */
    if(go_left){
      left->add_father(id);
    } else {
      right->add_father(id);
    } 
  }
  return true;
}

// [[Rcpp::export]]
NumericVector MP_build2(CharacterMatrix seqs, Function print){
  int nc = seqs.ncol();
  
  //构建初始树
  Ind *first = new Ind(-1, NULL, seqs);
  first->childs[0] = new Ind(0, first, seqs);
  first->childs[1] = new Ind(1, first, seqs);
  
  int mark = 0;
  Ind *new_top = NULL;
  //逐步添加
  for(int i = 2;i < nc;i++){
    //将加入序列置顶
    new_top = first->add_father(i);
    mark = cal_tree(new_top);
    if(place_the_seq(first, first, i, mark, print)){
      first->delete_father();
    } else{
      first = new_top;
    }
  }
  /*
  for(int i=0;i<seqs.nrow();i++){
    print('d');
    print(first->childs[0]->branch_base(i));
    print(first->childs[0]->branch_long);
    print(first->childs[1]->branch_base(i));
    print(first->childs[1]->branch_long);
    print(first->branch_base(i));
    print(first->branch_long);
    print('d');
  }*/
  
  /*
  int cnt = 10;
  int best_mark = cal_tree(first);
  int mark2 = 0;
  int best_mark2 = cal_tree2(first);
  double Tp = 100;
  print(best_mark);
  while(Tp > 1){
    cnt = 10;
    while(cnt > 0){
      first->update_num(0);
      int num1 = first->childs[0]->child_num+1;
      int num2 = first->childs[1]->child_num+1;
      Ind *l1[num1];
      Ind *l2[num2];
      first->childs[0]->front_ergodic(0, l1);
      first->childs[1]->front_ergodic(0, l2);
      Ind *ind1 = l1[rand() % num1];
      Ind *ind2 = l2[rand() % num2];
      if(!change_node(ind1, ind2)){
        continue;
      }
      mark = cal_tree(first);
      mark2 = cal_tree2(first);
      int deta =  best_mark - mark;
      //(get_random() < exp((double)deta/Tp) && Tp > 10)
      if(deta > 0 || (get_random() < exp((double)deta/Tp) && Tp > 10)){
        best_mark2 = best_mark2 < mark2 ? best_mark2 : mark2;
        best_mark = best_mark < mark ? best_mark : mark;
        cnt = deta > 0 ? 10 : cnt;
        if(deta <= 0){
          //print(best_mark2 - mark2);
        }
      } else {
        change_node(ind1, ind2);
      }
      cnt--; 
    }
    Tp *= 0.99;
  }
  print(best_mark);
  */
  
  
  /*
  int cnt = 1;
  int best_mark2 = cal_tree2(first), mark2 = 0;
  while(cnt > 0){
    first->update_num(0);
    first->right_lerp();
    first->left_mesh(first, print);
    mark2 = cal_tree2(first);
    if(mark2 >= best_mark2){
      break;
      best_mark2 = mark2;
    }
    cnt++;
  }*/
  
  //first->update_num(0);
  //first->right_lerp();
  //first->left_mesh(first, print);
  //first->right_lerp();
  //first->left_mesh(first, print);
  /*print(cal_tree(first));
  first->right_lerp();
  first->left_mesh(first, print);
  print(cal_tree(first));*/
  
  int apoint = nc + nc - 1;
  Ind *back_list[apoint];
  first->front_ergodic(0, back_list);
  NumericVector ans(apoint);
  for(int i =0;i<apoint;i++){
    ans[i] = back_list[i]->leaf_id;
  }
  delete first;
  return ans;
}

class Node{
public:
  int leaf_id;
  NumericVector dis;
  Node *childs[2];
  double mesh_range(Node *another);
  int front_ergodic(Node *back_list[], int st);
  Node(int id, NumericVector tdis, Node *child1, Node*child2);
  ~Node();
};

Node::Node(int id, NumericVector tdis, Node *child1, Node*child2){
  leaf_id = id;
  dis = tdis;
  childs[0] = child1;
  childs[1] = child2;
}

Node::~Node(){
  for(int i=0;i<2;i++){
    if(childs[i]!=NULL){
      delete childs[i];
    }
  }
}

double Node::mesh_range(Node *another){
  if(another==this) return 0;
  if(leaf_id != -1){
    if(another->leaf_id != -1){
      return dis[another->leaf_id];
    } else {
      return another->mesh_range(this);
    }
  }
  return (childs[0]->mesh_range(another) + childs[1]->mesh_range(another) - childs[0]->mesh_range(childs[1])) / 2;
}


int Node::front_ergodic(Node *back_list[], int st){
  if(leaf_id==-1){
    st = childs[0]->front_ergodic(back_list, st);
    st = childs[1]->front_ergodic(back_list, st);
  }
  back_list[st] = this; 
  st++;
  return st;
}

Node *mesh_node(Node *left, Node *right){
  Node *new_node = new Node(-1, NULL, left, right);
  return new_node;
}

double cal_node(int num, Node *wait[]){
  double S0 = 0;
  double still = 0;
  for(int j=0;j<num;j++){
    if(wait[j]!=NULL) still+=1;
    for(int k=j;k<num;k++){
      if(wait[j]==NULL||wait[k]==NULL||k==j) continue;
      S0 += wait[j]->mesh_range(wait[k]);
    }
  }
  if(still==1) return 0;
  return S0/(still-1);
}

// [[Rcpp::export]]
NumericVector NJ_build(NumericMatrix distance, Function print){
  int num = distance.ncol();
  
  Node *wait[num];
  for(int i=0;i<num;i++){
    wait[i] = new Node(i, distance(_,i), NULL, NULL);
  }
  Node *final_node = NULL;
  NumericVector r(num);
  double minS = -1;
  double tmpS = -1;
  int pl = 0;
  int pr = 0;
  int still = num;
  for(int i=1;i<num;i++){
    minS = -1;
    for(int u=0;u<num;u++){
      r[u] = 0;
      if(wait[u]==NULL) continue;
      for(int v=0;v<num;v++){
        if(wait[v]==NULL) continue;
        r[u] += wait[u]->mesh_range(wait[v]);
      }
    }
    
    for(int j=0;j<(num-1);j++){
      for(int k=(j+1);k<num;k++){
        if(wait[j]==NULL||wait[k]==NULL) continue;
        if(still==2){
          pl = j;
          pr = k;
          break;
        }
        /*
        Node *new_node = mesh_node(wait[j], wait[k]);
        wait[j] = new_node;
        wait[k] = NULL;
        tmpS = cal_node(num, wait);
        if(tmpS < minS || minS == -1){
          minS = tmpS;
          pl = j;
          pr = k;
        }
        wait[j] = new_node->childs[0];
        wait[k] = new_node->childs[1];
        new_node->childs[0] = NULL;
        new_node->childs[1] = NULL;
        delete new_node;
         */
        tmpS = wait[j]->mesh_range(wait[k]) / ((r[j]+r[j])/(still-2));
        if(tmpS < minS || minS == -1){
          minS = tmpS;
          pl = j;
          pr = k;
        }
      }
    }
    
    final_node = mesh_node(wait[pl], wait[pr]);
    wait[pl] = final_node;
    wait[pr] = NULL;
    still--;
  }
  
  int anode = num+num-1;
  NumericVector ans(anode);
  Node *back_list[anode];
  final_node->front_ergodic(back_list,0);
  for(int i=0;i<anode;i++){
    ans[i] = back_list[i]->leaf_id;
  }
  delete final_node;
  return ans;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
