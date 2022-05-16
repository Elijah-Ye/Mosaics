/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */

    if(first[curDim] > second[curDim])
    {
      return false;
    }

    if(first[curDim] < second[curDim])
    {
      return true;
    }

    return first < second;

}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */

    double dist1 = 0.0; //currentBest
    double dist2 = 0.0; //potential

    for(int i = 0; i < Dim; ++i){
      dist1 += (currentBest[i] - target[i])*(currentBest[i] - target[i]);
      dist2 += (potential[i] - target[i])*(potential[i] - target[i]);
    }

    if(dist2 < dist1)
    {
      return true;
    }
    
    if(dist2 > dist1)
    {
      return false;
    }
    
    return potential < currentBest;


}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
  /**
   * @todo Implement this function!
   */
    vector<Point<Dim>> Lists = newPoints;
    root = KDTreeHelper(Lists, 0, Lists.size() - 1, 0);

}

template <int Dim>
typename KDTree<Dim>::KDTreeNode * KDTree<Dim>::KDTreeHelper(vector<Point<Dim>>&Lists, unsigned left, unsigned right, int dim)
{
  size = 0;
  if (Lists.empty() || left > right || left<0 || right<0 || left >= Lists.size() || right>=Lists.size()) {
    return NULL;
  }
  unsigned mid = (left + right)/2;
  KDTreeNode * curRoot = new KDTreeNode(select(Lists, left, right, mid, dim % Dim));
  dim++; 
  size++;
  curRoot->left = KDTreeHelper(Lists, left, mid - 1, dim % Dim);
  curRoot->right = KDTreeHelper(Lists, mid + 1, right, dim % Dim);
  return curRoot;
}



//select
template <int Dim>
Point<Dim> KDTree<Dim>::select(vector<Point<Dim>>&Lists, unsigned left, unsigned right, unsigned k, int dim)
{
  if(left == right)
  {
    return Lists[left];
  }
  unsigned pivotIndex = partition(Lists, left, right, k, dim);

  if(k == pivotIndex)
  {
    return Lists[k];
  }

  else if(k < pivotIndex)
  { 
    return select(Lists, left, pivotIndex-1, k, dim);
  }

  else
  {
    return select(Lists, pivotIndex+1, right, k, dim);
  }

}

//partition
template <int Dim>
unsigned KDTree<Dim>::partition(vector<Point<Dim>>&Lists, unsigned left, unsigned right, unsigned pivotIndex, int dim)
{
  Point<Dim> pivotValue = Lists[pivotIndex];
  Point<Dim> temp;
  
  //swap list[pivotIndex] and list[right]
  temp = Lists[pivotIndex];
  Lists[pivotIndex] = Lists[right]; //move the pivot to the end
  Lists[right] = temp;

  unsigned storeIndex = left;

  for(unsigned i = left; i < right; i++) //I did not include the last element. Do we include the last element?
  { 
    if(smallerDimVal(Lists[i], pivotValue, dim)){

      //swap list[storeIndex] and list[i]
      temp = Lists[storeIndex];
      Lists[storeIndex] = Lists[i];
      Lists[i] = temp;

      storeIndex++;
    }
  }

  //swap list[right] and list[storeIndex]
  temp = Lists[right];
  Lists[right] = Lists[storeIndex];
  Lists[storeIndex] = temp;

  return storeIndex;
}



template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  size = other.size;
  copy(this->root, other.root);

}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode * root1, KDTreeNode * root2){
  if(root2 == nullptr)
  {
    root1 == nullptr;
  }
  else
  {
    root1 = new KDTreeNode();
    root1->point = root2->point;
    copy(root1->left, root2->left);
    copy(root1->right, root2->right);
  }
}


template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  if(this->root != rhs.root)
  {
    clear(this->root);
    copy(this->root, rhs.root);
  }

  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
  clear(this->root);
}

template <int Dim>
void KDTree<Dim>::clear(KDTreeNode *& root){
  if(root != nullptr)
  {
    clear(root->left);
    clear(root->right);
    delete root;
  }
}


template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    
    return findNearestNeighbor(query, 0, root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query, int dimension, KDTreeNode * curRoot) const{
  
  if(curRoot->left == nullptr && curRoot->right == nullptr)
  {
    return curRoot->point;
  }

  bool left = false; //use to check if we already check left subtree. 

  Point<Dim> nearest = curRoot->point;
  Point<Dim> tempNearest = nearest;

  if (smallerDimVal(query, curRoot->point, dimension) && curRoot->left != nullptr)
  {
    tempNearest = findNearestNeighbor(query, (dimension+1)%Dim, curRoot->left);
    left = true;
  }
  else if((!smallerDimVal(query, curRoot->point, dimension)) && curRoot->right != nullptr)
  {
    tempNearest = findNearestNeighbor(query, (dimension+1)%Dim, curRoot->right);
  }

  if(shouldReplace(query, nearest, tempNearest))
  {
    nearest = tempNearest;
  }  

  double radius = 0.0;
  double splitDist = 0.0;

  for(int i = 0; i < Dim; i++){
    radius += (query[i] - tempNearest[i])*(query[i] - tempNearest[i]);
  }

  splitDist = (curRoot->point[dimension] - query[dimension])*(curRoot->point[dimension] - query[dimension]);

  if(radius >= splitDist)
  {
    if(left && curRoot->right != nullptr)
    {
      tempNearest = findNearestNeighbor(query, (dimension+1)%Dim, curRoot->right);
    }
    else if(!left && curRoot->left != nullptr)
    {
      tempNearest = findNearestNeighbor(query, (dimension+1)%Dim, curRoot->left);
    }
    if(shouldReplace(query, nearest, tempNearest))
    {
      nearest = tempNearest;
    }
  }

  return nearest;
}



