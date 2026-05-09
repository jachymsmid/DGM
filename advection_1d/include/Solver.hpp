#pragma once

/**
 * @class Solver
 * @brief Manages all the parts of the DGM solver
 */
template < class Model >
class Solver
{
public:

   Solver() = default;

   void init(int argc, char* argv[]) {}

   Real compute_time_step() {}

   void initMesh() {}

   void initMeshData() {}

   void integrate() {}

   void makeSnapshot() {}

   void readUserConfig() {}

   void writeLog() {}

private:


};
