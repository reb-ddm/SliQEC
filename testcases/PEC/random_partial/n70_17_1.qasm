OPENQASM 2.0;
include "qelib1.inc";
qreg q[70];
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
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
cx q[33], q[27];
s q[24];
h q[8];
cx q[8], q[69];
tdg q[69];
cx q[47], q[69];
t q[69];
cx q[8], q[69];
tdg q[69];
cx q[47], q[69];
t q[69];
cx q[47], q[8];
tdg q[8];
cx q[47], q[8];
t q[47];
t q[8];
h q[8];
t q[56];
s q[63];
cx q[9], q[22];
t q[47];
cx q[12], q[2];
h q[18];
h q[37];
cx q[37], q[53];
tdg q[53];
cx q[24], q[53];
t q[53];
cx q[37], q[53];
tdg q[53];
cx q[24], q[53];
t q[53];
cx q[24], q[37];
tdg q[37];
cx q[24], q[37];
t q[24];
t q[37];
h q[37];
s q[42];
h q[38];
cx q[38], q[15];
tdg q[15];
cx q[49], q[15];
t q[15];
cx q[38], q[15];
tdg q[15];
cx q[49], q[15];
t q[15];
cx q[49], q[38];
tdg q[38];
cx q[49], q[38];
t q[49];
t q[38];
h q[38];
s q[22];
cx q[30], q[29];
t q[13];
t q[31];
h q[27];
cx q[27], q[50];
tdg q[50];
cx q[60], q[50];
t q[50];
cx q[27], q[50];
tdg q[50];
cx q[60], q[50];
t q[50];
cx q[60], q[27];
tdg q[27];
cx q[60], q[27];
t q[60];
t q[27];
h q[27];
h q[63];
cx q[63], q[15];
tdg q[15];
cx q[21], q[15];
t q[15];
cx q[63], q[15];
tdg q[15];
cx q[21], q[15];
t q[15];
cx q[21], q[63];
tdg q[63];
cx q[21], q[63];
t q[21];
t q[63];
h q[63];
h q[1];
cx q[1], q[38];
tdg q[38];
cx q[5], q[38];
t q[38];
cx q[1], q[38];
tdg q[38];
cx q[5], q[38];
t q[38];
cx q[5], q[1];
tdg q[1];
cx q[5], q[1];
t q[5];
t q[1];
h q[1];
h q[29];
cx q[29], q[8];
tdg q[8];
cx q[11], q[8];
t q[8];
cx q[29], q[8];
tdg q[8];
cx q[11], q[8];
t q[8];
cx q[11], q[29];
tdg q[29];
cx q[11], q[29];
t q[11];
t q[29];
h q[29];
h q[13];
h q[41];
h q[46];
cx q[46], q[49];
tdg q[49];
cx q[25], q[49];
t q[49];
cx q[46], q[49];
tdg q[49];
cx q[25], q[49];
t q[49];
cx q[25], q[46];
tdg q[46];
cx q[25], q[46];
t q[25];
t q[46];
h q[46];
cx q[60], q[37];
cx q[39], q[23];
cx q[45], q[25];
h q[55];
cx q[55], q[16];
tdg q[16];
cx q[15], q[16];
t q[16];
cx q[55], q[16];
tdg q[16];
cx q[15], q[16];
t q[16];
cx q[15], q[55];
tdg q[55];
cx q[15], q[55];
t q[15];
t q[55];
h q[55];
h q[42];
h q[68];
cx q[68], q[43];
tdg q[43];
cx q[27], q[43];
t q[43];
cx q[68], q[43];
tdg q[43];
cx q[27], q[43];
t q[43];
cx q[27], q[68];
tdg q[68];
cx q[27], q[68];
t q[27];
t q[68];
h q[68];
t q[18];
cx q[19], q[41];
h q[65];
h q[7];
cx q[7], q[39];
tdg q[39];
cx q[47], q[39];
t q[39];
cx q[7], q[39];
tdg q[39];
cx q[47], q[39];
t q[39];
cx q[47], q[7];
tdg q[7];
cx q[47], q[7];
t q[47];
t q[7];
h q[7];
h q[60];
t q[41];
h q[29];
cx q[29], q[37];
tdg q[37];
cx q[20], q[37];
t q[37];
cx q[29], q[37];
tdg q[37];
cx q[20], q[37];
t q[37];
cx q[20], q[29];
tdg q[29];
cx q[20], q[29];
t q[20];
t q[29];
h q[29];
h q[54];
cx q[54], q[21];
tdg q[21];
cx q[10], q[21];
t q[21];
cx q[54], q[21];
tdg q[21];
cx q[10], q[21];
t q[21];
cx q[10], q[54];
tdg q[54];
cx q[10], q[54];
t q[10];
t q[54];
h q[54];
s q[4];
cx q[57], q[42];
t q[42];
h q[11];
cx q[11], q[25];
tdg q[25];
cx q[62], q[25];
t q[25];
cx q[11], q[25];
tdg q[25];
cx q[62], q[25];
t q[25];
cx q[62], q[11];
tdg q[11];
cx q[62], q[11];
t q[62];
t q[11];
h q[11];
h q[50];
t q[48];
s q[25];
h q[43];
s q[51];
h q[37];
cx q[37], q[54];
tdg q[54];
cx q[56], q[54];
t q[54];
cx q[37], q[54];
tdg q[54];
cx q[56], q[54];
t q[54];
cx q[56], q[37];
tdg q[37];
cx q[56], q[37];
t q[56];
t q[37];
h q[37];
h q[52];
h q[59];
s q[8];
h q[31];
s q[47];
cx q[24], q[66];
cx q[31], q[40];
cx q[10], q[23];
h q[6];
cx q[6], q[26];
tdg q[26];
cx q[9], q[26];
t q[26];
cx q[6], q[26];
tdg q[26];
cx q[9], q[26];
t q[26];
cx q[9], q[6];
tdg q[6];
cx q[9], q[6];
t q[9];
t q[6];
h q[6];
h q[38];
cx q[38], q[44];
tdg q[44];
cx q[11], q[44];
t q[44];
cx q[38], q[44];
tdg q[44];
cx q[11], q[44];
t q[44];
cx q[11], q[38];
tdg q[38];
cx q[11], q[38];
t q[11];
t q[38];
h q[38];
h q[29];
t q[12];
h q[33];
h q[12];
cx q[12], q[50];
tdg q[50];
cx q[41], q[50];
t q[50];
cx q[12], q[50];
tdg q[50];
cx q[41], q[50];
t q[50];
cx q[41], q[12];
tdg q[12];
cx q[41], q[12];
t q[41];
t q[12];
h q[12];
t q[8];
h q[28];
cx q[28], q[12];
tdg q[12];
cx q[14], q[12];
t q[12];
cx q[28], q[12];
tdg q[12];
cx q[14], q[12];
t q[12];
cx q[14], q[28];
tdg q[28];
cx q[14], q[28];
t q[14];
t q[28];
h q[28];
s q[38];
s q[11];
s q[42];
h q[67];
cx q[53], q[40];
cx q[64], q[45];
cx q[16], q[37];
s q[35];
h q[48];
cx q[48], q[40];
tdg q[40];
cx q[66], q[40];
t q[40];
cx q[48], q[40];
tdg q[40];
cx q[66], q[40];
t q[40];
cx q[66], q[48];
tdg q[48];
cx q[66], q[48];
t q[66];
t q[48];
h q[48];
s q[38];
h q[4];
cx q[4], q[22];
tdg q[22];
cx q[47], q[22];
t q[22];
cx q[4], q[22];
tdg q[22];
cx q[47], q[22];
t q[22];
cx q[47], q[4];
tdg q[4];
cx q[47], q[4];
t q[47];
t q[4];
h q[4];
h q[45];
cx q[45], q[40];
tdg q[40];
cx q[17], q[40];
t q[40];
cx q[45], q[40];
tdg q[40];
cx q[17], q[40];
t q[40];
cx q[17], q[45];
tdg q[45];
cx q[17], q[45];
t q[17];
t q[45];
h q[45];
h q[68];
h q[54];
cx q[54], q[48];
tdg q[48];
cx q[3], q[48];
t q[48];
cx q[54], q[48];
tdg q[48];
cx q[3], q[48];
t q[48];
cx q[3], q[54];
tdg q[54];
cx q[3], q[54];
t q[3];
t q[54];
h q[54];
cx q[42], q[13];
s q[54];
h q[52];
cx q[13], q[26];
t q[11];
h q[1];
h q[29];
cx q[29], q[8];
tdg q[8];
cx q[18], q[8];
t q[8];
cx q[29], q[8];
tdg q[8];
cx q[18], q[8];
t q[8];
cx q[18], q[29];
tdg q[29];
cx q[18], q[29];
t q[18];
t q[29];
h q[29];
t q[29];
h q[12];
cx q[12], q[44];
tdg q[44];
cx q[26], q[44];
t q[44];
cx q[12], q[44];
tdg q[44];
cx q[26], q[44];
t q[44];
cx q[26], q[12];
tdg q[12];
cx q[26], q[12];
t q[26];
t q[12];
h q[12];
cx q[2], q[43];
h q[0];
s q[5];
cx q[33], q[55];
t q[14];
s q[49];
h q[10];
t q[39];
h q[46];
h q[47];
cx q[47], q[6];
tdg q[6];
cx q[56], q[6];
t q[6];
cx q[47], q[6];
tdg q[6];
cx q[56], q[6];
t q[6];
cx q[56], q[47];
tdg q[47];
cx q[56], q[47];
t q[56];
t q[47];
h q[47];
cx q[45], q[59];
h q[55];
h q[39];
cx q[39], q[4];
tdg q[4];
cx q[32], q[4];
t q[4];
cx q[39], q[4];
tdg q[4];
cx q[32], q[4];
t q[4];
cx q[32], q[39];
tdg q[39];
cx q[32], q[39];
t q[32];
t q[39];
h q[39];
h q[28];
cx q[20], q[58];
h q[17];
cx q[11], q[53];
h q[10];
cx q[10], q[17];
tdg q[17];
cx q[43], q[17];
t q[17];
cx q[10], q[17];
tdg q[17];
cx q[43], q[17];
t q[17];
cx q[43], q[10];
tdg q[10];
cx q[43], q[10];
t q[43];
t q[10];
h q[10];
cx q[23], q[12];
s q[14];
s q[0];
h q[21];
cx q[21], q[20];
tdg q[20];
cx q[17], q[20];
t q[20];
cx q[21], q[20];
tdg q[20];
cx q[17], q[20];
t q[20];
cx q[17], q[21];
tdg q[21];
cx q[17], q[21];
t q[17];
t q[21];
h q[21];
s q[6];
h q[13];
h q[63];
cx q[63], q[12];
tdg q[12];
cx q[50], q[12];
t q[12];
cx q[63], q[12];
tdg q[12];
cx q[50], q[12];
t q[12];
cx q[50], q[63];
tdg q[63];
cx q[50], q[63];
t q[50];
t q[63];
h q[63];
t q[68];
h q[57];
cx q[57], q[4];
tdg q[4];
cx q[25], q[4];
t q[4];
cx q[57], q[4];
tdg q[4];
cx q[25], q[4];
t q[4];
cx q[25], q[57];
tdg q[57];
cx q[25], q[57];
t q[25];
t q[57];
h q[57];
cx q[8], q[17];
h q[4];
h q[19];
cx q[19], q[27];
tdg q[27];
cx q[4], q[27];
t q[27];
cx q[19], q[27];
tdg q[27];
cx q[4], q[27];
t q[27];
cx q[4], q[19];
tdg q[19];
cx q[4], q[19];
t q[4];
t q[19];
h q[19];
s q[53];
h q[59];
cx q[59], q[29];
tdg q[29];
cx q[12], q[29];
t q[29];
cx q[59], q[29];
tdg q[29];
cx q[12], q[29];
t q[29];
cx q[12], q[59];
tdg q[59];
cx q[12], q[59];
t q[12];
t q[59];
h q[59];
t q[30];
t q[28];
h q[31];
h q[7];
cx q[7], q[54];
tdg q[54];
cx q[29], q[54];
t q[54];
cx q[7], q[54];
tdg q[54];
cx q[29], q[54];
t q[54];
cx q[29], q[7];
tdg q[7];
cx q[29], q[7];
t q[29];
t q[7];
h q[7];
s q[28];
t q[55];
h q[28];
h q[4];
cx q[4], q[45];
tdg q[45];
cx q[65], q[45];
t q[45];
cx q[4], q[45];
tdg q[45];
cx q[65], q[45];
t q[45];
cx q[65], q[4];
tdg q[4];
cx q[65], q[4];
t q[65];
t q[4];
h q[4];
h q[48];
cx q[48], q[57];
tdg q[57];
cx q[14], q[57];
t q[57];
cx q[48], q[57];
tdg q[57];
cx q[14], q[57];
t q[57];
cx q[14], q[48];
tdg q[48];
cx q[14], q[48];
t q[14];
t q[48];
h q[48];
t q[52];
h q[47];
h q[8];
h q[11];
h q[0];
t q[37];
t q[7];
h q[20];
cx q[20], q[19];
tdg q[19];
cx q[7], q[19];
t q[19];
cx q[20], q[19];
tdg q[19];
cx q[7], q[19];
t q[19];
cx q[7], q[20];
tdg q[20];
cx q[7], q[20];
t q[7];
t q[20];
h q[20];
h q[45];
h q[39];
t q[16];
t q[39];
cx q[69], q[39];
s q[23];
cx q[45], q[8];
h q[65];
h q[45];
cx q[45], q[33];
tdg q[33];
cx q[46], q[33];
t q[33];
cx q[45], q[33];
tdg q[33];
cx q[46], q[33];
t q[33];
cx q[46], q[45];
tdg q[45];
cx q[46], q[45];
t q[46];
t q[45];
h q[45];
h q[11];
cx q[11], q[36];
tdg q[36];
cx q[15], q[36];
t q[36];
cx q[11], q[36];
tdg q[36];
cx q[15], q[36];
t q[36];
cx q[15], q[11];
tdg q[11];
cx q[15], q[11];
t q[15];
t q[11];
h q[11];
s q[7];
s q[43];
h q[20];
cx q[20], q[62];
tdg q[62];
cx q[49], q[62];
t q[62];
cx q[20], q[62];
tdg q[62];
cx q[49], q[62];
t q[62];
cx q[49], q[20];
tdg q[20];
cx q[49], q[20];
t q[49];
t q[20];
h q[20];
s q[8];
h q[69];
s q[39];
cx q[56], q[58];
cx q[46], q[10];
t q[51];
h q[13];
cx q[13], q[55];
tdg q[55];
cx q[15], q[55];
t q[55];
cx q[13], q[55];
tdg q[55];
cx q[15], q[55];
t q[55];
cx q[15], q[13];
tdg q[13];
cx q[15], q[13];
t q[15];
t q[13];
h q[13];
s q[24];
h q[13];
cx q[13], q[28];
tdg q[28];
cx q[26], q[28];
t q[28];
cx q[13], q[28];
tdg q[28];
cx q[26], q[28];
t q[28];
cx q[26], q[13];
tdg q[13];
cx q[26], q[13];
t q[26];
t q[13];
h q[13];
h q[15];
cx q[15], q[57];
t q[13];
s q[43];
s q[9];
h q[40];
cx q[8], q[47];
cx q[65], q[15];
h q[41];
cx q[41], q[18];
tdg q[18];
cx q[68], q[18];
t q[18];
cx q[41], q[18];
tdg q[18];
cx q[68], q[18];
t q[18];
cx q[68], q[41];
tdg q[41];
cx q[68], q[41];
t q[68];
t q[41];
h q[41];
cx q[20], q[27];
h q[60];
h q[22];
cx q[22], q[11];
tdg q[11];
cx q[24], q[11];
t q[11];
cx q[22], q[11];
tdg q[11];
cx q[24], q[11];
t q[11];
cx q[24], q[22];
tdg q[22];
cx q[24], q[22];
t q[24];
t q[22];
h q[22];
h q[21];
cx q[21], q[65];
tdg q[65];
cx q[55], q[65];
t q[65];
cx q[21], q[65];
tdg q[65];
cx q[55], q[65];
t q[65];
cx q[55], q[21];
tdg q[21];
cx q[55], q[21];
t q[55];
t q[21];
h q[21];
s q[15];
h q[33];
t q[56];
cx q[63], q[55];
h q[58];
h q[1];
cx q[1], q[64];
tdg q[64];
cx q[4], q[64];
t q[64];
cx q[1], q[64];
tdg q[64];
cx q[4], q[64];
t q[64];
cx q[4], q[1];
tdg q[1];
cx q[4], q[1];
t q[4];
t q[1];
h q[1];
t q[46];
t q[25];
h q[12];
t q[28];
cx q[56], q[58];
t q[39];
h q[53];
cx q[53], q[40];
tdg q[40];
cx q[6], q[40];
t q[40];
cx q[53], q[40];
tdg q[40];
cx q[6], q[40];
t q[40];
cx q[6], q[53];
tdg q[53];
cx q[6], q[53];
t q[6];
t q[53];
h q[53];
cx q[62], q[60];
h q[56];
t q[14];
h q[66];
h q[62];
s q[64];
cx q[16], q[52];
cx q[24], q[7];
s q[6];
t q[26];
h q[38];
cx q[38], q[17];
tdg q[17];
cx q[16], q[17];
t q[17];
cx q[38], q[17];
tdg q[17];
cx q[16], q[17];
t q[17];
cx q[16], q[38];
tdg q[38];
cx q[16], q[38];
t q[16];
t q[38];
h q[38];
cx q[54], q[31];
s q[2];
cx q[43], q[36];
s q[4];
t q[17];
s q[35];
h q[13];
cx q[11], q[9];
s q[29];
t q[45];
h q[60];
s q[21];
s q[60];
h q[46];
s q[14];
h q[8];
cx q[8], q[25];
tdg q[25];
cx q[22], q[25];
t q[25];
cx q[8], q[25];
tdg q[25];
cx q[22], q[25];
t q[25];
cx q[22], q[8];
tdg q[8];
cx q[22], q[8];
t q[22];
t q[8];
h q[8];
z q[0];
t q[1];
cx q[1], q[2];
x q[3];
cx q[4], q[5];
t q[5];
cx q[4], q[5];
x q[7];
y q[12];
x q[14];
y q[19];
cx q[23], q[22];
x q[22];
s q[22];
cx q[23], q[22];
s q[25];
cx q[24], q[25];
x q[26];
cx q[28], q[27];
s q[30];
z q[31];
y q[32];
h q[38];
cx q[38], q[53];
tdg q[53];
cx q[58], q[53];
t q[53];
cx q[38], q[53];
tdg q[53];
cx q[58], q[53];
t q[53];
cx q[58], q[38];
tdg q[38];
cx q[58], q[38];
t q[58];
t q[38];
h q[38];
h q[58];
cx q[58], q[68];
tdg q[68];
cx q[50], q[68];
t q[68];
cx q[58], q[68];
tdg q[68];
cx q[50], q[68];
t q[68];
cx q[50], q[58];
tdg q[58];
cx q[50], q[58];
t q[50];
t q[58];
h q[58];
t q[50];
s q[68];
s q[36];
h q[40];
h q[47];
h q[46];
cx q[46], q[35];
tdg q[35];
cx q[68], q[35];
t q[35];
cx q[46], q[35];
tdg q[35];
cx q[68], q[35];
t q[35];
cx q[68], q[46];
tdg q[46];
cx q[68], q[46];
t q[68];
t q[46];
h q[46];
h q[66];
h q[62];
cx q[62], q[49];
tdg q[49];
cx q[48], q[49];
t q[49];
cx q[62], q[49];
tdg q[49];
cx q[48], q[49];
t q[49];
cx q[48], q[62];
tdg q[62];
cx q[48], q[62];
t q[48];
t q[62];
h q[62];
h q[38];
h q[64];
h q[45];
cx q[45], q[64];
tdg q[64];
cx q[37], q[64];
t q[64];
cx q[45], q[64];
tdg q[64];
cx q[37], q[64];
t q[64];
cx q[37], q[45];
tdg q[45];
cx q[37], q[45];
t q[37];
t q[45];
h q[45];
t q[47];
s q[46];
h q[42];
cx q[42], q[61];
tdg q[61];
cx q[54], q[61];
t q[61];
cx q[42], q[61];
tdg q[61];
cx q[54], q[61];
t q[61];
cx q[54], q[42];
tdg q[42];
cx q[54], q[42];
t q[54];
t q[42];
h q[42];
h q[63];
h q[63];
cx q[63], q[35];
tdg q[35];
cx q[43], q[35];
t q[35];
cx q[63], q[35];
tdg q[35];
cx q[43], q[35];
t q[35];
cx q[43], q[63];
tdg q[63];
cx q[43], q[63];
t q[43];
t q[63];
h q[63];
s q[46];
h q[66];
h q[59];
h q[62];
s q[67];
h q[53];
cx q[53], q[58];
tdg q[58];
cx q[62], q[58];
t q[58];
cx q[53], q[58];
tdg q[58];
cx q[62], q[58];
t q[58];
cx q[62], q[53];
tdg q[53];
cx q[62], q[53];
t q[62];
t q[53];
h q[53];
t q[45];
t q[36];
cx q[46], q[64];
h q[35];
s q[50];
h q[42];
s q[66];
cx q[67], q[62];
h q[35];
s q[49];
cx q[41], q[57];
