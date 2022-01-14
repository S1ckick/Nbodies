// #include <fstream>
// #include <iostream>

// #include "pointmasses.h"
// int main() {
//   std::ifstream infile("../NBodies/points.txt");
//   long double x, y, z, vx, vy, vz, m;
//   std::string name;

//   long double masses[16];
//   long double dx[96];
//   long double dy[96];
//   long double dz[96];
//   long double dist[96];
//   long double dist2[96];
//   long double dist3[96];
//   long double temp_acc[96];
//   long double arr[96];
//   long double f[96];
//   int i = 0;
//   while (infile >> name >> m >> x >> y >> z >> vx >> vy >> vz) {
//     masses[i] = m;
//     arr[i*6] = x;
//     arr[i * 6 + 1] = y;
//     arr[i * 6 + 2] = z;
//     arr[i * 6 + 3] = vx;
//     arr[i * 6 + 4] = vy;
//     arr[i * 6 + 5] = vz;

//     std::cout << "parsed " << name << " : " << x << " " << y << " " << z << " "
//               << vx << " " << vy << " " << vz << " " << m << std::endl;
//     i++;
//   }

//   PointMasses pm;
//   pm.n_objects = 16;

//   pm.masses = masses;
//   pm.dx = dx;
//   pm.dy = dy;
//   pm.dz = dz;
//   pm.dist = dist;
//   pm.dist2 = dist2;
//   pm.dist3 = dist3;
//   pm.temp_acc = temp_acc;

//   pointmassesCalculateXdot(arr, f, &pm);

//   for (int j = 0; j < 16; j++) {
//     printf("x: %f, y: %f, z: %f, vx: %f, vy: %f, vz: %f \n", f[j * 6],
//            f[j * 6 + 1], f[j * 6 + 2], f[j * 6 + 3], f[j * 6 + 4],
//            f[j * 6 + 5]);
//   }

//   return 0;
// }