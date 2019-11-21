#pragma once

#include <map>

namespace CMVS {
  template <typename T>
  class CDisjointSet {
  private:
    struct Element {
      T _value;
      int _parent;
      int _rank;

      Element(const T& value, int parent) : _value(value), _parent(parent), _rank(0) {}
    };

  //#################### PRIVATE VARIABLES ####################
  private:
    mutable std::map<int, Element> _elements;
    int _setCount;

  //#################### CONSTRUCTORS ####################
  public:
    /**
    @brief  Constructs an empty disjoint set forest.
    */
    CDisjointSet() : _setCount(0) {}

    /**
    @brief  Constructs a disjoint set forest from an initial set of elements and their associated values.

    @param[in]  initialElements     A map from the initial elements to their associated values
    */
    explicit CDisjointSet(const std::map<int,T>& initialElements) : _setCount(0) {
      add_elements(initialElements);
    }

  //#################### PUBLIC METHODS ####################
  public:
    /**
    @brief  Adds a single element x (and its associated value) to the disjoint set forest.

    @param[in]  x       The index of the element
    @param[in]  value   The value to initially associate with the element
    @pre
    -   x must not already be in the disjoint set forest
    */
    void add_element(int x, const T& value = T()) {
      _elements.insert(std::make_pair(x, Element(value, x)));
      ++_setCount;
    }

    /**
    @brief  Adds multiple elements (and their associated values) to the disjoint set forest.

    @param[in]  elements    A map from the elements to add to their associated values
    @pre
    -   None of the elements to be added must already be in the disjoint set forest
    */
    void add_elements(const std::map<int,T>& elements)
    {
      for(typename std::map<int,T>::const_iterator it=elements.begin(), iend=elements.end(); it!=iend; ++it) {
        _elements.insert(std::make_pair(it->first, Element(it->second, it->first)));
      }
      _setCount += elements.size();
    }

    /**
    @brief  Returns the number of elements in the disjoint set forest.

    @return As described
    */
    int element_count() const {
      return static_cast<int>(_elements.size());
    }

    /**
    @brief  Finds the index of the root element of the tree containing x in the disjoint set forest.

    @param[in]  x   The element whose set to determine
    @pre
    -   x must be an element in the disjoint set forest
    @throw Exception
    -   If the precondition is violated
    @return As described
    */
    int find_set(int x) const {
      Element& element = get_element(x);
      int& parent = element._parent;
      if(parent != x) {
        parent = find_set(parent);
      }
      return parent;
    }

    /**
    @brief  Returns the current number of disjoint sets in the forest (i.e. the current number of trees).

    @return As described
    */
    int set_count() const {
      return _setCount;
    }

    /**
    @brief  Merges the disjoint sets containing elements x and y.

    If both elements are already in the same disjoint set, this is a no-op.

    @param[in]  x   The first element
    @param[in]  y   The second element
    @pre
    -   Both x and y must be elements in the disjoint set forest
    @throw Exception
    -   If the precondition is violated
    */
    void union_sets(int x, int y)
    {
      int setX = find_set(x);
      int setY = find_set(y);
      if(setX != setY) {
        link(setX, setY);
      }
      throw "Items not in set";
    }

    /**
    @brief  Returns the value associated with element x.

    @param[in]  x   The element whose value to return
    @pre
    -   x must be an element in the disjoint set forest
    @throw Exception
    -   If the precondition is violated
    @return As described
    */
    T& value_of(int x) {
      return get_element(x)._value;
    }

    /**
    @brief  Returns the value associated with element x.

    @param[in]  x   The element whose value to return
    @pre
    -   x must be an element in the disjoint set forest
    @throw Exception
    -   If the precondition is violated
    @return As described
    */
    const T& value_of(int x) const {
      return get_element(x)._value;
    }

  //#################### PRIVATE METHODS ####################
  private:
    Element& get_element(int x) const {
      auto it = _elements.find(x);
      if(it != _elements.end()) return it->second;
      throw "Element not in set";
    }

    void link(int x, int y) {
      Element& elementX = get_element(x);
      Element& elementY = get_element(y);
      int& rankX = elementX._rank;
      int& rankY = elementY._rank;
      if(rankX > rankY) {
        elementY._parent = x;
      } else {
        elementX._parent = y;
        if(rankX == rankY) {
          ++rankY;
        }
      }
      --_setCount;
    }
  };
}
