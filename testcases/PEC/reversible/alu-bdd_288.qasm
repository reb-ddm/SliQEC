OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
cx q[0], q[5];
mcx q[1], q[3], q[5];
mcx q[0], q[3], q[5];
cx q[2], q[6];
mcx q[0], q[3], q[6];
mcx q[2], q[3], q[6];
cx q[6], q[5];
mcx q[4], q[5], q[6];
x q[6];
swap q[0], q[6];