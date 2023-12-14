#include <iostream>
#include <cstring>

char* in() {
  char* str = new char;
  int symbol;
  int len = 0;
  int str_capacity = 2;
  symbol = getchar();
  while (!isspace(symbol)) {
    if (len + 1 >= str_capacity) {
      str_capacity *= 2;
      char* str1 = new char[str_capacity];
      memcpy(str1, str, len);
      delete str;
      str = str1;
    }
    str[len] = static_cast<char>(symbol);
    ++len;
    symbol = getchar();
  }
  str[len] = '\0';
  return str;
}

void push(char* elem, char**& stack, int& size, int& capacity) {
  if (size + 1 >= capacity) {
    capacity *= 2;
    char** new_stack = new char* [capacity];
    for (int i = 0; i <= size; ++i) {
      new_stack[i] = new char;
    }
    memcpy(new_stack, stack, size);
    delete[] stack;
    stack = new_stack;
  }
  stack[size] = elem;
  ++size;
}

const char* pop(char**& stack, int& size) {
  if (size == 0) {
    return static_cast<const char*>("error");
  }
  char** new_stack = new char* [size - 1];
  --size;
  for (int i = 0; i < size; ++i) {
    new_stack[i] = stack[i];
  }
  char* ans = stack[size];
  delete[] stack;
  stack = new_stack;
  return ans;
}

const char* back(char** stack, int size) {
  if (size == 0) {
    return static_cast<const char*>("error");
  }
  return stack[size - 1];
}

void clear(char**& stack, int& size) {
  char** new_stack = new char* [0];
  delete[] stack;
  stack = new_stack;
  size = 0;
}


int main() {
  int size = 0;
  int capacity = 2;
  char** stack = new char* [0];
  while(true){
    char* query = in();
    if (strcmp(query, "push") == 0) {
      char* str = in();
      push(str, stack, size, capacity);
      delete str;
      std::cout << "ok" << '\n';
      continue;
    }
    if (strcmp(query,"pop") == 0) {
      std::cout << pop(stack, size) << '\n';
      continue;
    }
    if (strcmp(query,"back") == 0) {
      std::cout << back(stack, size) << '\n';
      continue;
    }
    if (strcmp(query,"size") == 0) {
      std::cout << size << '\n';
      continue;
    }
    if(strcmp(query,"clear") == 0) {
      clear(stack, size);
      std::cout << "ok" << '\n';
      continue;
    }
    if(strcmp(query, "exit") == 0) {
      for (int i = 0; i < size; ++i) {
        delete stack[i];
      }
      delete[] stack;
      std::cout << "bye";
      break;
    }
    delete query;
  }
  return 0;
}
