#include <gtest/gtest.h>
#include "../headers/Mesh.hpp"

static constexpr double TOL = 1e-12;

TEST(MeshTest, UniformMeshElementCount)
{
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 10);
    EXPECT_EQ(mesh.numElements(), 10);
}

TEST(MeshTest, UniformMeshFaceCount)
{
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 10);
    EXPECT_EQ(mesh.numFaces(), 11);
}

TEST(MeshTest, UniformMeshVertexCoordinates)
{
    int K = 5;
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
    for (int k = 0; k < K; ++k) {
        EXPECT_NEAR(mesh.leftVertex(k),  double(k)   / K, TOL);
        EXPECT_NEAR(mesh.rightVertex(k), double(k+1) / K, TOL);
    }
}

TEST(MeshTest, UniformMeshElementSize)
{
    int K = 4;
    auto mesh = DG::Mesh<double>::uniform(0.0, 2.0, K);
    for (int k = 0; k < K; ++k)
        EXPECT_NEAR(mesh.elementSize(k), 0.5, TOL);
}

TEST(MeshTest, JacobianIsHalfElementSize)
{
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 4);
    for (int k = 0; k < 4; ++k)
        EXPECT_NEAR(mesh.jacobian(k), mesh.elementSize(k) * 0.5, TOL);
}

TEST(MeshTest, BoundaryFaceDetection)
{
    int K = 5;
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
    EXPECT_TRUE(mesh.isBoundaryFace(0));
    EXPECT_TRUE(mesh.isBoundaryFace(K));
    for (int f = 1; f < K; ++f)
        EXPECT_FALSE(mesh.isBoundaryFace(f));
}

TEST(MeshTest, FaceNeighbours)
{
    int K = 5;
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, K);
    // Internal face f=2: left=1, right=2
    EXPECT_EQ(mesh.leftCellOfFace(2),  1);
    EXPECT_EQ(mesh.rightCellOfFace(2), 2);
    // Left boundary
    EXPECT_EQ(mesh.leftCellOfFace(0),  DG::Mesh<double>::BOUNDARY_FACE);
    // Right boundary
    EXPECT_EQ(mesh.rightCellOfFace(K), DG::Mesh<double>::BOUNDARY_FACE);
}

TEST(MeshTest, Normals)
{
    auto mesh = DG::Mesh<double>::uniform(0.0, 1.0, 4);
    EXPECT_DOUBLE_EQ(mesh.leftNormal(),  -1.0);
    EXPECT_DOUBLE_EQ(mesh.rightNormal(), +1.0);
}
