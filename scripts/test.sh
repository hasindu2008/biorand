#!/bin/bash

BINARY="biorand"
TEST_PATH="test/olp/"
	./$BINARY olp $TEST_PATH/t0.in $TEST_PATH/t0.out 15 2 5 10 2>   $TEST_PATH/err.txt && sdiff -b $TEST_PATH/t0.exp $TEST_PATH/t0.out
	./$BINARY olp $TEST_PATH/t1.in $TEST_PATH/t1.out 15 2 5 10 2>	 $TEST_PATH/err.txt && sdiff -b $TEST_PATH/t1.exp $TEST_PATH/t1.out
	./$BINARY olp $TEST_PATH/t3.in $TEST_PATH/t3.out 15 4 5 10 2>   $TEST_PATH/err.txt && sdiff -b $TEST_PATH/t3.exp $TEST_PATH/t3.out
	./$BINARY olp $TEST_PATH/t2.in $TEST_PATH/t2.out 30 2 15 25 2>  $TEST_PATH/err.txt && sdiff -b $TEST_PATH/t2.exp $TEST_PATH/t2.out
	./$BINARY olp $TEST_PATH/t4.in $TEST_PATH/t4.out 20 2 5 10 2>   $TEST_PATH/err.txt && sdiff -b $TEST_PATH/t4.exp $TEST_PATH/t4.out
	./$BINARY olp $TEST_PATH/t5.in $TEST_PATH/t5.out 125 3 15 50 2> $TEST_PATH/err.txt && sdiff -b $TEST_PATH/t5.exp $TEST_PATH/t5.out