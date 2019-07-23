#include <math.h>  /* sin */
#include <stdio.h> /* printf */
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
  double theta1, theta2, theta3;
  theta1 = 2 * M_PI;  // first rotation (yaw)
  theta2 = 2 * M_PI;  // second rotation (pitch)
  theta3 = 2 * M_PI;  // third rotation (roll)

  std::array<std::array<double, 3>, 3> rotation_3D{
      {{cos(theta1) * cos(theta2), cos(theta2) * sin(theta1), -sin(theta2)},
       {cos(theta1) * sin(theta2) * sin(theta3) - cos(theta3) * sin(theta1),
        cos(theta1) * cos(theta3) + sin(theta1) * sin(theta2) * sin(theta3),
        cos(theta2) * sin(theta3)},
       {cos(theta1) * cos(theta3) * sin(theta2) + sin(theta1) * sin(theta3),
        cos(theta3) * sin(theta1) * sin(theta2) - cos(theta1) * sin(theta3),
        cos(theta2) * cos(theta3)}}};

  std::array<std::array<double, 3>, 3> matrix_3x3{
      {{0.0, 0.0, 3.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};

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
  std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>
      matrix_new2;

  for (size_t ib = 0; ib < 3; ++ib) {
    for (size_t jb = 0; jb < 3; ++jb) {
      for (size_t kb = 0; kb < 3; ++kb) {
        for (size_t lb = 0; lb < 3; ++lb) {
          matrix_new2[ib][jb][kb][lb] = 0.0;
        }
      }
    }
  }

  std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> matrix_old;

  for (size_t ib = 0; ib < 3; ++ib) {
    for (size_t jb = 0; jb < 3; ++jb) {
      for (size_t kb = 0; kb < 3; ++kb) {
        for (size_t lb = 0; lb < 3; ++lb) {
          matrix_old[ib][jb][kb][lb] = rand() % 100;
        }
      }
    }
  }

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t ib = 0; ib < 3; ++ib) {
        for (size_t jb = 0; jb < 3; ++jb) {
          for (size_t k = 0; k < 3; ++k) {
            for (size_t kb = 0; kb < 3; ++kb) {
              for (size_t l = 0; l < 3; ++l) {
                for (size_t lb = 0; lb < 3; ++lb) {
                  matrix_new2[ib][jb][kb][lb] +=
                      rotation_3D[ib][i] * rotation_3D[jb][j] *
                      rotation_3D[kb][k] * rotation_3D[lb][l] *
                      matrix_old[i][j][k][l];
                }
              }
            }
          }
        }
      }
    }
  }

  print_tensor(matrix_old);
  std::cout << "\n\n";
  print_tensor(matrix_new2);
  return 0;
}

// clang-format -i -style=Google /Users/geoffrey/Desktop/Dot.cpp
// compile: clang++ -std=c++17 -o Dot Dot.cpp
