#include "../include/polyakov_loops.h"

#include <map>
#include <vector>

std::map<double, double>
polyakov_average_directions(const std::vector<double> &correlators, int D_max) {
  std::map<double, double> result;
  std::map<double, int> corr_num;

  double distance;
  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;
  int Dx_min;

  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;
      for (int Dx = Dx_min; Dx <= D_max; Dx++) {
        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
        if (!(Dy != 0 && Dz != 0 && Dx < 0)) {
          if (distance <= D_max && !(Dx == 0 && Dy == 0 && Dz == 0)) {
            result[distance] +=
                correlators[(Dz + D_max) * result_size1 +
                            (Dy + D_max) * result_size + Dx + D_max];
            corr_num[distance]++;
          }
        }
      }
    }
  }
  for (auto it = result.begin(); it != result.end(); ++it) {
    it->second /= corr_num[it->first];
  }

  return result;
}