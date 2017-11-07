#include "rng.h"

std::random_device rd;
std::mt19937 Gen(rd());
