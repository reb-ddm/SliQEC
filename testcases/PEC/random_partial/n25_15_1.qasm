OPENQASM 2.0;
include "qelib1.inc";
qreg q[25];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
s q[16];
cx q[3], q[17];
s q[12];
cx q[10], q[23];
h q[2];
h q[22];
h q[9];
cx q[9], q[17];
tdg q[17];
cx q[22], q[17];
t q[17];
cx q[9], q[17];
tdg q[17];
cx q[22], q[17];
t q[17];
cx q[22], q[9];
tdg q[9];
cx q[22], q[9];
t q[22];
t q[9];
h q[9];
cx q[5], q[3];
t q[19];
h q[16];
h q[23];
cx q[23], q[19];
tdg q[19];
cx q[4], q[19];
t q[19];
cx q[23], q[19];
tdg q[19];
cx q[4], q[19];
t q[19];
cx q[4], q[23];
tdg q[23];
cx q[4], q[23];
t q[4];
t q[23];
h q[23];
h q[0];
t q[3];
t q[21];
h q[21];
s q[11];
h q[6];
cx q[6], q[24];
tdg q[24];
cx q[20], q[24];
t q[24];
cx q[6], q[24];
tdg q[24];
cx q[20], q[24];
t q[24];
cx q[20], q[6];
tdg q[6];
cx q[20], q[6];
t q[20];
t q[6];
h q[6];
cx q[13], q[19];
h q[3];
h q[1];
s q[22];
s q[13];
t q[23];
cx q[22], q[5];
cx q[5], q[9];
s q[13];
s q[1];
cx q[23], q[10];
s q[13];
h q[15];
h q[2];
s q[11];
h q[0];
cx q[0], q[13];
tdg q[13];
cx q[12], q[13];
t q[13];
cx q[0], q[13];
tdg q[13];
cx q[12], q[13];
t q[13];
cx q[12], q[0];
tdg q[0];
cx q[12], q[0];
t q[12];
t q[0];
h q[0];
h q[12];
cx q[12], q[14];
tdg q[14];
cx q[5], q[14];
t q[14];
cx q[12], q[14];
tdg q[14];
cx q[5], q[14];
t q[14];
cx q[5], q[12];
tdg q[12];
cx q[5], q[12];
t q[5];
t q[12];
h q[12];
h q[2];
h q[18];
cx q[18], q[8];
tdg q[8];
cx q[9], q[8];
t q[8];
cx q[18], q[8];
tdg q[8];
cx q[9], q[8];
t q[8];
cx q[9], q[18];
tdg q[18];
cx q[9], q[18];
t q[9];
t q[18];
h q[18];
t q[11];
s q[13];
t q[18];
h q[9];
t q[20];
t q[7];
h q[13];
cx q[13], q[11];
tdg q[11];
cx q[2], q[11];
t q[11];
cx q[13], q[11];
tdg q[11];
cx q[2], q[11];
t q[11];
cx q[2], q[13];
tdg q[13];
cx q[2], q[13];
t q[2];
t q[13];
h q[13];
h q[11];
cx q[11], q[20];
tdg q[20];
cx q[21], q[20];
t q[20];
cx q[11], q[20];
tdg q[20];
cx q[21], q[20];
t q[20];
cx q[21], q[11];
tdg q[11];
cx q[21], q[11];
t q[21];
t q[11];
h q[11];
s q[4];
h q[21];
s q[3];
cx q[15], q[4];
cx q[7], q[10];
s q[18];
h q[21];
cx q[21], q[14];
tdg q[14];
cx q[16], q[14];
t q[14];
cx q[21], q[14];
tdg q[14];
cx q[16], q[14];
t q[14];
cx q[16], q[21];
tdg q[21];
cx q[16], q[21];
t q[16];
t q[21];
h q[21];
t q[1];
cx q[16], q[14];
cx q[20], q[9];
h q[6];
s q[16];
cx q[14], q[4];
s q[20];
h q[23];
cx q[23], q[18];
tdg q[18];
cx q[20], q[18];
t q[18];
cx q[23], q[18];
tdg q[18];
cx q[20], q[18];
t q[18];
cx q[20], q[23];
tdg q[23];
cx q[20], q[23];
t q[20];
t q[23];
h q[23];
h q[21];
s q[23];
h q[3];
cx q[3], q[15];
tdg q[15];
cx q[1], q[15];
t q[15];
cx q[3], q[15];
tdg q[15];
cx q[1], q[15];
t q[15];
cx q[1], q[3];
tdg q[3];
cx q[1], q[3];
t q[1];
t q[3];
h q[3];
t q[2];
s q[24];
t q[20];
s q[20];
h q[15];
h q[22];
cx q[22], q[1];
tdg q[1];
cx q[13], q[1];
t q[1];
cx q[22], q[1];
tdg q[1];
cx q[13], q[1];
t q[1];
cx q[13], q[22];
tdg q[22];
cx q[13], q[22];
t q[13];
t q[22];
h q[22];
h q[4];
cx q[4], q[5];
tdg q[5];
cx q[11], q[5];
t q[5];
cx q[4], q[5];
tdg q[5];
cx q[11], q[5];
t q[5];
cx q[11], q[4];
tdg q[4];
cx q[11], q[4];
t q[11];
t q[4];
h q[4];
cx q[18], q[2];
h q[3];
cx q[3], q[8];
tdg q[8];
cx q[13], q[8];
t q[8];
cx q[3], q[8];
tdg q[8];
cx q[13], q[8];
t q[8];
cx q[13], q[3];
tdg q[3];
cx q[13], q[3];
t q[13];
t q[3];
h q[3];
t q[2];
s q[12];
cx q[9], q[24];
s q[9];
y q[0];
t q[1];
cx q[2], q[1];
s q[1];
cx q[4], q[3];
s q[3];
cx q[5], q[6];
s q[10];
h q[12];
s q[23];
cx q[16], q[22];
s q[16];
h q[22];
cx q[22], q[18];
tdg q[18];
cx q[14], q[18];
t q[18];
cx q[22], q[18];
tdg q[18];
cx q[14], q[18];
t q[18];
cx q[14], q[22];
tdg q[22];
cx q[14], q[22];
t q[14];
t q[22];
h q[22];
t q[15];
t q[24];
cx q[16], q[20];
t q[22];
cx q[15], q[18];
s q[21];
s q[19];
t q[20];
