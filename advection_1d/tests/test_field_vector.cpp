/**
 * @file test_field_vector.cpp
 * @brief Unit tests for FieldVector.
 */

#include <gtest/gtest.h>
#include "FieldVector.hpp"

// ── Default constructor ────────────────────────────────────────────────────────
TEST(FieldVectorTest, DefaultConstructorSizes)
{
    DG::FieldVector<double> fv;
    EXPECT_EQ(fv.numElements(), 0);
    EXPECT_EQ(fv.numDOF(), 0);
    EXPECT_EQ(fv.totalSize(), 0);
}

// ── Parameterised constructor ──────────────────────────────────────────────────
TEST(FieldVectorTest, ParameterizedConstructorSizes)
{
    DG::FieldVector<double> fv(3, 4);
    EXPECT_EQ(fv.numElements(), 3);
    EXPECT_EQ(fv.numDOF(), 4);
    EXPECT_EQ(fv.totalSize(), 12);
}

TEST(FieldVectorTest, SingleElementSingleDOF)
{
    DG::FieldVector<double> fv(1, 1);
    EXPECT_EQ(fv.numElements(), 1);
    EXPECT_EQ(fv.numDOF(), 1);
    EXPECT_EQ(fv.totalSize(), 1);
}

// ── Zero initialisation ───────────────────────────────────────────────────────
TEST(FieldVectorTest, ZeroInitialization)
{
    DG::FieldVector<double> fv(3, 4);
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < 4; ++i)
            EXPECT_EQ(fv.elementPtr(k)[i], 0.0)
                << "Non-zero at element " << k << " dof " << i;
}

// ── elementPtr read/write ──────────────────────────────────────────────────────
TEST(FieldVectorTest, ElementPtrReadWrite)
{
    DG::FieldVector<double> fv(2, 3);
    fv.elementPtr(0)[0] = 1.0;
    fv.elementPtr(0)[1] = 2.0;
    fv.elementPtr(0)[2] = 3.0;
    fv.elementPtr(1)[0] = 4.0;
    fv.elementPtr(1)[1] = 5.0;
    fv.elementPtr(1)[2] = 6.0;

    EXPECT_EQ(fv.elementPtr(0)[0], 1.0);
    EXPECT_EQ(fv.elementPtr(0)[1], 2.0);
    EXPECT_EQ(fv.elementPtr(0)[2], 3.0);
    EXPECT_EQ(fv.elementPtr(1)[0], 4.0);
    EXPECT_EQ(fv.elementPtr(1)[1], 5.0);
    EXPECT_EQ(fv.elementPtr(1)[2], 6.0);
}

TEST(FieldVectorTest, ElementPtrContiguousLayout)
{
    // Element k starts at index k * Np in the underlying storage.
    DG::FieldVector<double> fv(3, 4);
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < 4; ++i)
            fv.elementPtr(k)[i] = static_cast<double>(k * 4 + i);

    const double* base = fv.data().getData();
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < 4; ++i)
            EXPECT_EQ(base[k * 4 + i], static_cast<double>(k * 4 + i));
}

// ── Const elementPtr ──────────────────────────────────────────────────────────
TEST(FieldVectorTest, ConstElementPtrReadable)
{
    DG::FieldVector<double> fv(2, 2);
    fv.elementPtr(0)[0] = 7.0;
    fv.elementPtr(1)[1] = 9.0;

    const DG::FieldVector<double>& cfv = fv;
    EXPECT_EQ(cfv.elementPtr(0)[0], 7.0);
    EXPECT_EQ(cfv.elementPtr(1)[1], 9.0);
}

// ── data() accessor ───────────────────────────────────────────────────────────
TEST(FieldVectorTest, DataAccessorSize)
{
    DG::FieldVector<double> fv(3, 5);
    EXPECT_EQ(fv.data().getSize(), 15);
}

TEST(FieldVectorTest, ConstDataAccessorSize)
{
    const DG::FieldVector<double> fv(4, 3);
    EXPECT_EQ(fv.data().getSize(), 12);
}

// ── copyFrom ──────────────────────────────────────────────────────────────────
TEST(FieldVectorTest, CopyFromValues)
{
    DG::FieldVector<double> src(2, 2);
    src.elementPtr(0)[0] = 1.0;
    src.elementPtr(0)[1] = 2.0;
    src.elementPtr(1)[0] = 3.0;
    src.elementPtr(1)[1] = 4.0;

    DG::FieldVector<double> dst(2, 2);
    dst.copyFrom(src);

    EXPECT_EQ(dst.elementPtr(0)[0], 1.0);
    EXPECT_EQ(dst.elementPtr(0)[1], 2.0);
    EXPECT_EQ(dst.elementPtr(1)[0], 3.0);
    EXPECT_EQ(dst.elementPtr(1)[1], 4.0);
}

TEST(FieldVectorTest, CopyFromIsDeep)
{
    // Modifying the source after copy must not affect the destination.
    DG::FieldVector<double> src(1, 2);
    src.elementPtr(0)[0] = 5.0;
    src.elementPtr(0)[1] = 6.0;

    DG::FieldVector<double> dst(1, 2);
    dst.copyFrom(src);

    src.elementPtr(0)[0] = 99.0;
    EXPECT_EQ(dst.elementPtr(0)[0], 5.0);
    EXPECT_EQ(dst.elementPtr(0)[1], 6.0);
}
