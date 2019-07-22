#include <array>
#include <iostream>

template <size_t Dim>
void print_array(std::array<double, Dim> vector) {
  std::cout << "[";
  for (size_t i = 0; i < vector.size(); ++i) {
    // std::cout << vector[i]; OK, but no bounds checking
    std::cout << vector.at(i);  // check i in bounds of vector
    if (i < vector.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]";
}
template <size_t Dim>
void print_tensor(std::array<double, Dim> tensor) {
  std::cout << "[";
  print_array(tensor);
  std::cout << "]\n";
}

template <size_t Dim>
void print_tensor(std::array<std::array<double, Dim>, Dim> tensor) {
  std::cout << "[";
  for (size_t i = 0; i < Dim; ++i) {
    print_array(tensor[i]);
    if (i < Dim - 1) {
      std::cout << ",\n";
    } else {
      std::cout << "]\n";
    }
  }
}

template <size_t Dim>
void print_tensor(
    std::array<std::array<std::array<double, Dim>, Dim>, Dim> tensor) {
  std::cout << "[";
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      print_array(tensor[i][j]);

      if (i < Dim - 1) {
        std::cout << ",\n";
      } else {
        std::cout << "]\n";
      }
    }
  }
}
template <size_t Dim>
void print_tensor(
    std::array<std::array<std::array<std::array<double, Dim>, Dim>, Dim>, Dim>
        tensor) {
  std::cout << "[";
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      for (size_t k = 0; k < Dim; ++k) {
        print_array(tensor[i][j][k]);
        if (i < Dim - 1) {
          std::cout << ",\n";
        } else {
          std::cout << "]\n";
        }
      }
    }
  }
}

int main() {
  std::array<double, 4> vector_a;
  vector_a[0] = 11111.0;
  vector_a[1] = 2.0;
  vector_a[2] = 1.0;
  vector_a[3] = -1.0;

  std::array<double, 3> vector_b{{2.0 / 3.0, 444.0, 4444.0}};

  std::array<std::array<double, 3>, 3> rotation_3D{
      {{0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}};

  std::array<std::array<double, 3>, 3> matrix_3x3{
      {{1.0, -1.0, 3.0}, {1.0, 2.0, 1.0}, {1.0, 1.0, 1.0}}};

  std::array<double, 3> vector_c{{0.0, 0.0, 0.0}};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      vector_c[i] += rotation_3D[i][j] * vector_b[j];
    }
  }

  std::array<std::array<double, 3>, 3> matrix_new;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t ib = 0; ib < 3; ++ib) {
        for (size_t jb = 0; jb < 3; ++jb) {
          matrix_new[ib][jb] +=
              rotation_3D[ib][i] * rotation_3D[jb][j] * matrix_3x3[i][j];
        }
      }
    }
  }

  std::array<std::array<std::array<double, 3>, 3>, 3> matrix_new1;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t ib = 0; ib < 3; ++ib) {
        for (size_t jb = 0; jb < 3; ++jb) {
          for (size_t k = 0; k < 3; ++k) {
            for (size_t kb = 0; kb < 3; ++kb) {
              matrix_new1[ib][jb][kb] += rotation_3D[ib][i] *
                                         rotation_3D[jb][j] *
                                         rotation_3D[kb][k] * matrix_3x3[i][j];
            }
          }
        }
      }
    }
  }

  std::array<double, 4> vector_nest;
  vector_a[0] = 0.0;
  vector_a[1] = 0.0;
  vector_a[2] = 0.0;
  vector_a[3] = 0.0;

  print_tensor(rotation_3D);
  print_tensor(vector_a);
  print_tensor(vector_b);
  print_tensor(vector_c);
  print_tensor(matrix_new1);

  return 0;
}

// clang-format -i -style=Google /Users/geoffrey/Desktop/Dot.cpp
// compile: clang++ -std=c++17 -o Dot Dot.cpp
