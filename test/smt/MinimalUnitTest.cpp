//
// Created by Tobi Popoola on 8/30/20.
//

#include <gtest/gtest.h>
#include <iegenlib.h>
class MinimalUnitTest : public::testing::Test{
protected:
    virtual void SetUp(){}
    virtual void TearDown(){}
};


TEST_F(MinimalUnitTest,CHECK_MINIM){
    // a < b + 1 equivalent to b - a >= 0

    // b - a
    Exp * e1 = new Exp();
    e1->addTerm(new VarTerm("b"));
    e1->addTerm(new VarTerm(-1,"a") );
    e1->setInequality();
    // b - a >= 0

    EXPECT_EQ("b - a >= 0", e1->toString());
}