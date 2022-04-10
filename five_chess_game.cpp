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

// [[Rcpp::export]]
int judge_game_over(IntegerMatrix board, int i,int j){
  int zi = board(i, j);
  if(zi){
    int cnt = -1;
    for(int it = i;it < 15;it++){
      if(board(it,j) != zi){
        break;
      }
      cnt++;
    }
    for(int it = i;it > 0;it--){
      if(board(it,j) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return zi;
    }
    
    cnt = -1;
    for(int jt = j;jt < 15;jt++){
      if(board(i,jt) != zi){
        break;
      }
      cnt++;
    }
    for(int jt = j;jt > 0;jt--){
      if(board(i,jt) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return zi;
    }
    
    cnt = -1;
    int im,jm;
    int mm = i < j ? i : j;
    for(int m = 0;m <= mm;m++){
      im = i - m;
      jm = j - m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    mm = 15 - mm;
    for(int m = 0;m < mm;m++){
      im = i + m;
      jm = j + m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return zi;
    }
    
    cnt = -1;
    mm = i < (14 - j) ? i : (14 - j);
    for(int m = 0;m <= mm;m++){
      im = i - m;
      jm = j + m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    mm = j < (14 - i) ? j : (14 - i);
    for(int m = 0;m < mm;m++){
      im = i + m;
      jm = j - m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return zi;
    }
  }
  return 0;
}


int put_chess(IntegerMatrix board, int i, int j, bool isDark) {
  int zi = isDark ? 1 : 2;
  if(board(i,j) == 0){
    board(i,j) = zi;
    
    int cnt = -1;
    for(int it = i;it < 15;it++){
      if(board(it,j) != zi){
        break;
      }
      cnt++;
    }
    for(int it = i;it > 0;it--){
      if(board(it,j) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return 1;
    }
    
    cnt = -1;
    for(int jt = j;jt < 15;jt++){
      if(board(i,jt) != zi){
        break;
      }
      cnt++;
    }
    for(int jt = j;jt > 0;jt--){
      if(board(i,jt) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return 1;
    }
    
    cnt = -1;
    int im,jm;
    int mm = i < j ? i : j;
    for(int m = 0;m <= mm;m++){
      im = i - m;
      jm = j - m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    mm = 15 - mm;
    for(int m = 0;m < mm;m++){
      im = i + m;
      jm = j + m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return 1;
    }
    
    cnt = -1;
    mm = i < (14 - j) ? i : (14 - j);
    for(int m = 0;m <= mm;m++){
      im = i - m;
      jm = j + m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    mm = j < (14 - i) ? j : (14 - i);
    for(int m = 0;m < mm;m++){
      im = i + m;
      jm = j - m;
      if(board(im, jm) != zi){
        break;
      }
      cnt++;
    }
    if(cnt > 4){
      return 1;
    }
    
    return 0;
  } else {
    return 2;
  }
}

// [[Rcpp::export]]
int random_game(IntegerMatrix board, bool isDark, IntegerVector spaces) {
  int num = spaces.size();
  int ans, i;
  int now;
  now = isDark;
  for(i = 0;i < num;i++){
    ans = put_chess(board, spaces[i] % 15, floor(spaces[i] / 15), now);
    if(ans){
      break;
    }
    now = !now;
  }
  ans = now ? 1 : 2;
  return ans;
}


class Node
{
public :
  Node *father;
  int win;
  int all;
  int zi;
  int deep;
  IntegerMatrix board;
  IntegerVector choices;
  Node *childs[225];
  int far;
  double UCT(void);
  void back(bool dark_win);
  int simulation(Function getSpaces);
  void MCTS(Function getAfter, Function getChoices, Function getSpaces);
  Node();
  void setValue(Node *fa);
  ~Node();
};

Node::Node(){
  win = 0;
  all = 0;
  far = 0;
  zi = 1;
  deep = 0;
}

void Node::setValue(Node *fa){
  if(fa != NULL){
    father = fa;
    deep = father->deep + 1;
    zi = father->zi == 1 ? 2 : 1;
  }
}

double Node::UCT(void){
  double xc = (double)win / (double)all;
  xc = zi==1 ? xc : (1 - xc);
  xc = xc + 0.6 * sqrt((double)(log(father->all)/all));
  return xc;
}

void Node::back(bool dark_win){
  all++;
  win = dark_win ? win + 1 : win;
  if(deep > 0){
    father->back(dark_win);
  }
}

int Node::simulation(Function getSpaces){
  IntegerVector spaces = getSpaces(board, zi==1);
  if(spaces.size()){
    int ans = random_game(board, zi==1, spaces);
    return ans == 1;
  }
  return true;
}

void Node::MCTS(Function getAfter, Function getChoices, Function getSpaces){
  int stop = choices.size();
  if(far < stop){
    childs[far] = new Node();
    childs[far]->setValue(this);
    childs[far]->board = getAfter(board, choices[far], zi);
    childs[far]->choices = getChoices(childs[far]->board, zi!=1, deep);
    bool win = childs[far]->simulation(getSpaces);
    childs[far]->back(win);
    far++;
  } else {
    int _far = far - 1;
    double max = childs[_far]->UCT();
    int w_max = _far;
    double tmp;
    for(int i=0;i<_far;i++){
      if(childs[i]!=NULL){
        tmp = childs[i]->UCT();
        if(tmp > max){
          max = tmp;
          w_max = i;
        }
      }
    }
    childs[w_max]->MCTS(getAfter, getChoices, getSpaces);
  }
}

Node::~Node(void){
  for(int i=0;i<far;i++){
    if(childs[i] != NULL){
      delete childs[i];
    }
  }
}

// [[Rcpp::export]]
int MCTS(IntegerMatrix board, bool isDark, Function getAfter, Function getChoices, Function getSpaces, Function debug, int times = 100000){
  Node root;
  root.zi = isDark ? 1 : 2;
  root.board = board;
  root.choices = getChoices(board, isDark, 0);
  for(int i = 0; i<times; i++){
    root.MCTS(getAfter, getChoices, getSpaces);
  }
  int far = root.far - 1;
  int max = root.childs[far]->all;
  int w_max = far;
  int tmp;
  for(int i = 0;i<far;i++){
    tmp = root.childs[i]->all;
    debug(tmp);
    if(tmp > max){
      max = tmp;
      w_max = i;
    }
  }
  return root.choices[w_max];
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
