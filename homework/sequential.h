//
// Created by simon on 25/06/2022.
//

#ifndef HOMEWORK_SEQUENTIAL_H
#define HOMEWORK_SEQUENTIAL_H

pc::matrix<int> floydWarhsallSeq(pc::matrix<double> &dist);
pc::matrix<int> floydWarhsallSquared(pc::matrix<double> &dist, int s);
void floyd(pc::matrix<double> &dist, pc::matrix<int> &pred, int dx, int dy, int ax, int ay, int bx, int by, int s);

#endif //HOMEWORK_SEQUENTIAL_H
