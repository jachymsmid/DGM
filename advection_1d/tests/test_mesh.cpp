#include <gtest/gtest.h>
#include "../headers/Mesh.hpp"

// set tolerance
static constexpr double TOL = 1e-12;

// test uniform constructor: number of elements
TEST(UniformMeshTest, UniformMeshElementCount)
{
  auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 10);
  EXPECT_EQ(mesh.numElements(), 10);
}

// test uniform constructor: number of faces
TEST(UniformMeshTest, UniformMeshFaceCount)
{
  auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 10);
  EXPECT_EQ(mesh.numFaces(), 11);
}

// test uniform constructor: vertex coordinates
TEST(UniformMeshTest, UniformMeshVertexCoordinates)
{
  int K = 5;
  auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
  for (int k = 0; k < K; ++k)
  {
    EXPECT_NEAR(mesh.leftVertex(k), double(k) / K, TOL);
    EXPECT_NEAR(mesh.rightVertex(k), double(k+1) / K, TOL);
  }
}

// test uniform constructor: vertex coordinates of neighbouring elements
TEST(UniformMeshTest, UniformMeshVertexCoordinatesNeighbours)
{
  int K = 5;
  auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
  for (int k = 0; k < K; ++k)
  {
    EXPECT_EQ(mesh.leftVertex(k + 1), mesh.rightVertex(k));
  }
}

// test uniform constructor: element size
TEST(UniformMeshTest, UniformMeshElementSize)
{
  int K = 4;
  auto mesh = DG::Mesh<double>::uniform(0.0, 2.0, K);
  for (int k = 0; k < K; ++k)
    EXPECT_NEAR(mesh.elementSize(k), 2.0/double(K), TOL);
}

// test uniform constructor: jacobian
TEST(UniformMeshTest, Jacobian)
{
  int K = 4;
  auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
  for (int k = 0; k < K; ++k)
    EXPECT_NEAR(mesh.jacobian(k), mesh.elementSize(k) * 0.5, TOL);
}

// test uniform constructor: is face boundary
TEST(UniformMeshTest, BoundaryFaceDetection)
{
  int K = 5;
  auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
  EXPECT_TRUE(mesh.isBoundaryFace(0));
  EXPECT_TRUE(mesh.isBoundaryFace(K));
  for (int f = 1; f < K; ++f)
    EXPECT_FALSE(mesh.isBoundaryFace(f));
}

// test uniform constructor: neighbours of faces
TEST(UniformMeshTest, FaceNeighbours)
{
    int K = 5;
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);

    // Internal face f=2: left=1, right=2
    EXPECT_EQ(mesh.leftCellOfFace(2),  1);
    EXPECT_EQ(mesh.rightCellOfFace(2), 2);

    // Left boundary
    EXPECT_EQ(mesh.leftCellOfFace(0),  DG::Mesh<double>::BOUNDARY_FACE);
    // Right boundary
    //
    EXPECT_EQ(mesh.rightCellOfFace(K), DG::Mesh<double>::BOUNDARY_FACE);
}

// test uniform constructor: normals at faces
TEST(UniformMeshTest, Normals)
{
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 4);
    EXPECT_DOUBLE_EQ(mesh.leftNormal(),  -1.0);
    EXPECT_DOUBLE_EQ(mesh.rightNormal(), +1.0);
}
