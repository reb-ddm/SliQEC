OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
cx q[3], q[4];
x q[3];
cx q[1], q[2];
mcx q[3], q[2], q[1];
cx q[1], q[0];
x q[1];
mcx q[4], q[0], q[1];
swap q[0], q[1];