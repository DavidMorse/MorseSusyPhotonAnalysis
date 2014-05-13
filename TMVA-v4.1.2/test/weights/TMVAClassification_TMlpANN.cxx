#include "weights/TMVAClassification_TMlpANN.h"
#include <cmath>

double TMVAClassification_TMlpANN::Value(int index,double in0,double in1,double in2,double in3,double in4,double in5,double in6) {
   input0 = (in0 - 1407.68)/1007.06;
   input1 = (in1 - 2.44929)/0.790649;
   input2 = (in2 - 1121.9)/727.381;
   input3 = (in3 - 151.844)/246.38;
   input4 = (in4 - 229.95)/291.02;
   input5 = (in5 - 137.17)/231.155;
   input6 = (in6 - 1.62389)/0.903669;
   switch(index) {
     case 0:
         return neuron0x87e4e60();
     default:
         return 0.;
   }
}

double TMVAClassification_TMlpANN::Value(int index, double* input) {
   input0 = (input[0] - 1407.68)/1007.06;
   input1 = (input[1] - 2.44929)/0.790649;
   input2 = (input[2] - 1121.9)/727.381;
   input3 = (input[3] - 151.844)/246.38;
   input4 = (input[4] - 229.95)/291.02;
   input5 = (input[5] - 137.17)/231.155;
   input6 = (input[6] - 1.62389)/0.903669;
   switch(index) {
     case 0:
         return neuron0x87e4e60();
     default:
         return 0.;
   }
}

double TMVAClassification_TMlpANN::neuron0x81cf310() {
   return input0;
}

double TMVAClassification_TMlpANN::neuron0x6d57820() {
   return input1;
}

double TMVAClassification_TMlpANN::neuron0x6d57ad0() {
   return input2;
}

double TMVAClassification_TMlpANN::neuron0x8e3bc50() {
   return input3;
}

double TMVAClassification_TMlpANN::neuron0x8e3bf90() {
   return input4;
}

double TMVAClassification_TMlpANN::neuron0x8d11a60() {
   return input5;
}

double TMVAClassification_TMlpANN::neuron0x8d11da0() {
   return input6;
}

double TMVAClassification_TMlpANN::input0x8066360() {
   double input = -1.63683;
   input += synapse0x8e3bbe0();
   input += synapse0x8066610();
   input += synapse0x8066650();
   input += synapse0x8066690();
   input += synapse0x80666d0();
   input += synapse0x8066710();
   input += synapse0x8066750();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8066360() {
   double input = input0x8066360();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8066790() {
   double input = -0.474003;
   input += synapse0x6d57f40();
   input += synapse0x887f5e0();
   input += synapse0x887f620();
   input += synapse0x8001150();
   input += synapse0x8001190();
   input += synapse0x8e3d8f0();
   input += synapse0x8e3d930();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8066790() {
   double input = input0x8066790();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8e492f0() {
   double input = -0.143153;
   input += synapse0x8e49630();
   input += synapse0x8e49670();
   input += synapse0x8e496b0();
   input += synapse0x89e73b0();
   input += synapse0x89e73f0();
   input += synapse0x8e3c2d0();
   input += synapse0x8e3c310();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8e492f0() {
   double input = input0x8e492f0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8e49800() {
   double input = 0.392927;
   input += synapse0x8d12170();
   input += synapse0x8d121b0();
   input += synapse0x8e3c350();
   input += synapse0x8e3c390();
   input += synapse0x849fe20();
   input += synapse0x849fe60();
   input += synapse0x849fea0();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8e49800() {
   double input = input0x8e49800();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x849fee0() {
   double input = 0.491573;
   input += synapse0x84a0220();
   input += synapse0x84a0260();
   input += synapse0x84a02a0();
   input += synapse0x84a02e0();
   input += synapse0x84a0320();
   input += synapse0x81747b0();
   input += synapse0x7fee620();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x849fee0() {
   double input = input0x849fee0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8edb310() {
   double input = -0.556933;
   input += synapse0x87ed270();
   input += synapse0x8edb650();
   input += synapse0x8edb690();
   input += synapse0x8edb6d0();
   input += synapse0x8edb710();
   input += synapse0x8edb750();
   input += synapse0x8edb790();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8edb310() {
   double input = input0x8edb310();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8edb7d0() {
   double input = -0.152448;
   input += synapse0x8edbb10();
   input += synapse0x87f2750();
   input += synapse0x74c4d50();
   input += synapse0x8806e60();
   input += synapse0x8806ea0();
   input += synapse0x8806ee0();
   input += synapse0x84a0570();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8edb7d0() {
   double input = input0x8edb7d0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x87f2dd0() {
   double input = 0.359657;
   input += synapse0x87f3110();
   input += synapse0x87f3150();
   input += synapse0x87f3190();
   input += synapse0x87f31d0();
   input += synapse0x87f3210();
   input += synapse0x87f3250();
   input += synapse0x87f3290();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x87f2dd0() {
   double input = input0x87f2dd0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x87f32d0() {
   double input = 0.878107;
   input += synapse0x84a05b0();
   input += synapse0x84a05f0();
   input += synapse0x8e496f0();
   input += synapse0x8e49730();
   input += synapse0x8e49770();
   input += synapse0x8e497b0();
   input += synapse0x8e40440();
   input += synapse0x8e40480();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x87f32d0() {
   double input = input0x87f32d0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8e404c0() {
   double input = 0.150619;
   input += synapse0x7b9d640();
   input += synapse0x7b9d680();
   input += synapse0x84a0360();
   input += synapse0x84a03a0();
   input += synapse0x84a03e0();
   input += synapse0x84a0420();
   input += synapse0x84a0460();
   input += synapse0x84a04a0();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8e404c0() {
   double input = input0x8e404c0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8e40910() {
   double input = 0.140675;
   input += synapse0x8e40bc0();
   input += synapse0x8e40c00();
   input += synapse0x8e40c40();
   input += synapse0x82dcf00();
   input += synapse0x82dcf40();
   input += synapse0x82dcf80();
   input += synapse0x82dcfc0();
   input += synapse0x82dd000();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8e40910() {
   double input = input0x8e40910();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x82dd040() {
   double input = -0.193147;
   input += synapse0x8337660();
   input += synapse0x83376a0();
   input += synapse0x83376e0();
   input += synapse0x8337720();
   input += synapse0x8337760();
   input += synapse0x83377a0();
   input += synapse0x83377e0();
   input += synapse0x8337820();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x82dd040() {
   double input = input0x82dd040();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x8337860() {
   double input = 0.456901;
   input += synapse0x8337ba0();
   input += synapse0x8337be0();
   input += synapse0x8337c20();
   input += synapse0x8337c60();
   input += synapse0x8337ca0();
   input += synapse0x8337ce0();
   input += synapse0x8337d20();
   input += synapse0x8337d60();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x8337860() {
   double input = input0x8337860();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x715e200() {
   double input = -0.169695;
   input += synapse0x715e540();
   input += synapse0x715e580();
   input += synapse0x715e5c0();
   input += synapse0x715e600();
   input += synapse0x715e640();
   input += synapse0x715e680();
   input += synapse0x715e6c0();
   input += synapse0x715e700();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x715e200() {
   double input = input0x715e200();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x715e740() {
   double input = 0.667952;
   input += synapse0x82dd2f0();
   input += synapse0x8337da0();
   input += synapse0x8337de0();
   input += synapse0x87e4d20();
   input += synapse0x87e4d60();
   input += synapse0x87e4da0();
   input += synapse0x87e4de0();
   input += synapse0x87e4e20();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x715e740() {
   double input = input0x715e740();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double TMVAClassification_TMlpANN::input0x87e4e60() {
   double input = -0.720768;
   input += synapse0x87e51a0();
   input += synapse0x87e51e0();
   input += synapse0x87e5220();
   input += synapse0x87e5260();
   input += synapse0x87e52a0();
   input += synapse0x87e52e0();
   input += synapse0x87e5320();
   return input;
}

double TMVAClassification_TMlpANN::neuron0x87e4e60() {
   double input = input0x87e4e60();
   return (input * 1)+0;
}

double TMVAClassification_TMlpANN::synapse0x8e3bbe0() {
   return (neuron0x81cf310()*-2.55601);
}

double TMVAClassification_TMlpANN::synapse0x8066610() {
   return (neuron0x6d57820()*0.610492);
}

double TMVAClassification_TMlpANN::synapse0x8066650() {
   return (neuron0x6d57ad0()*-0.455809);
}

double TMVAClassification_TMlpANN::synapse0x8066690() {
   return (neuron0x8e3bc50()*-2.02909);
}

double TMVAClassification_TMlpANN::synapse0x80666d0() {
   return (neuron0x8e3bf90()*-4.43217);
}

double TMVAClassification_TMlpANN::synapse0x8066710() {
   return (neuron0x8d11a60()*-7.95347);
}

double TMVAClassification_TMlpANN::synapse0x8066750() {
   return (neuron0x8d11da0()*0.707843);
}

double TMVAClassification_TMlpANN::synapse0x6d57f40() {
   return (neuron0x81cf310()*-1.14371);
}

double TMVAClassification_TMlpANN::synapse0x887f5e0() {
   return (neuron0x6d57820()*0.697362);
}

double TMVAClassification_TMlpANN::synapse0x887f620() {
   return (neuron0x6d57ad0()*-0.834067);
}

double TMVAClassification_TMlpANN::synapse0x8001150() {
   return (neuron0x8e3bc50()*-0.394779);
}

double TMVAClassification_TMlpANN::synapse0x8001190() {
   return (neuron0x8e3bf90()*-0.298035);
}

double TMVAClassification_TMlpANN::synapse0x8e3d8f0() {
   return (neuron0x8d11a60()*-0.397864);
}

double TMVAClassification_TMlpANN::synapse0x8e3d930() {
   return (neuron0x8d11da0()*0.101202);
}

double TMVAClassification_TMlpANN::synapse0x8e49630() {
   return (neuron0x81cf310()*-1.54184);
}

double TMVAClassification_TMlpANN::synapse0x8e49670() {
   return (neuron0x6d57820()*0.0624082);
}

double TMVAClassification_TMlpANN::synapse0x8e496b0() {
   return (neuron0x6d57ad0()*-0.952498);
}

double TMVAClassification_TMlpANN::synapse0x89e73b0() {
   return (neuron0x8e3bc50()*-1.83782);
}

double TMVAClassification_TMlpANN::synapse0x89e73f0() {
   return (neuron0x8e3bf90()*-1.73769);
}

double TMVAClassification_TMlpANN::synapse0x8e3c2d0() {
   return (neuron0x8d11a60()*-2.303);
}

double TMVAClassification_TMlpANN::synapse0x8e3c310() {
   return (neuron0x8d11da0()*0.644465);
}

double TMVAClassification_TMlpANN::synapse0x8d12170() {
   return (neuron0x81cf310()*1.08712);
}

double TMVAClassification_TMlpANN::synapse0x8d121b0() {
   return (neuron0x6d57820()*-0.491231);
}

double TMVAClassification_TMlpANN::synapse0x8e3c350() {
   return (neuron0x6d57ad0()*0.814014);
}

double TMVAClassification_TMlpANN::synapse0x8e3c390() {
   return (neuron0x8e3bc50()*0.992675);
}

double TMVAClassification_TMlpANN::synapse0x849fe20() {
   return (neuron0x8e3bf90()*1.89906);
}

double TMVAClassification_TMlpANN::synapse0x849fe60() {
   return (neuron0x8d11a60()*4.27591);
}

double TMVAClassification_TMlpANN::synapse0x849fea0() {
   return (neuron0x8d11da0()*-0.719364);
}

double TMVAClassification_TMlpANN::synapse0x84a0220() {
   return (neuron0x81cf310()*-1.93102);
}

double TMVAClassification_TMlpANN::synapse0x84a0260() {
   return (neuron0x6d57820()*0.419557);
}

double TMVAClassification_TMlpANN::synapse0x84a02a0() {
   return (neuron0x6d57ad0()*-0.670559);
}

double TMVAClassification_TMlpANN::synapse0x84a02e0() {
   return (neuron0x8e3bc50()*-3.91936);
}

double TMVAClassification_TMlpANN::synapse0x84a0320() {
   return (neuron0x8e3bf90()*-4.35569);
}

double TMVAClassification_TMlpANN::synapse0x81747b0() {
   return (neuron0x8d11a60()*-2.77456);
}

double TMVAClassification_TMlpANN::synapse0x7fee620() {
   return (neuron0x8d11da0()*0.178017);
}

double TMVAClassification_TMlpANN::synapse0x87ed270() {
   return (neuron0x81cf310()*-0.210747);
}

double TMVAClassification_TMlpANN::synapse0x8edb650() {
   return (neuron0x6d57820()*0.193125);
}

double TMVAClassification_TMlpANN::synapse0x8edb690() {
   return (neuron0x6d57ad0()*-0.0164421);
}

double TMVAClassification_TMlpANN::synapse0x8edb6d0() {
   return (neuron0x8e3bc50()*0.0247438);
}

double TMVAClassification_TMlpANN::synapse0x8edb710() {
   return (neuron0x8e3bf90()*-0.0307957);
}

double TMVAClassification_TMlpANN::synapse0x8edb750() {
   return (neuron0x8d11a60()*-0.237853);
}

double TMVAClassification_TMlpANN::synapse0x8edb790() {
   return (neuron0x8d11da0()*0.0188079);
}

double TMVAClassification_TMlpANN::synapse0x8edbb10() {
   return (neuron0x81cf310()*1.75685);
}

double TMVAClassification_TMlpANN::synapse0x87f2750() {
   return (neuron0x6d57820()*-0.626355);
}

double TMVAClassification_TMlpANN::synapse0x74c4d50() {
   return (neuron0x6d57ad0()*0.998752);
}

double TMVAClassification_TMlpANN::synapse0x8806e60() {
   return (neuron0x8e3bc50()*1.38129);
}

double TMVAClassification_TMlpANN::synapse0x8806ea0() {
   return (neuron0x8e3bf90()*1.68387);
}

double TMVAClassification_TMlpANN::synapse0x8806ee0() {
   return (neuron0x8d11a60()*3.8312);
}

double TMVAClassification_TMlpANN::synapse0x84a0570() {
   return (neuron0x8d11da0()*-0.126857);
}

double TMVAClassification_TMlpANN::synapse0x87f3110() {
   return (neuron0x81cf310()*1.07321);
}

double TMVAClassification_TMlpANN::synapse0x87f3150() {
   return (neuron0x6d57820()*-0.491291);
}

double TMVAClassification_TMlpANN::synapse0x87f3190() {
   return (neuron0x6d57ad0()*0.703034);
}

double TMVAClassification_TMlpANN::synapse0x87f31d0() {
   return (neuron0x8e3bc50()*0.949614);
}

double TMVAClassification_TMlpANN::synapse0x87f3210() {
   return (neuron0x8e3bf90()*0.629716);
}

double TMVAClassification_TMlpANN::synapse0x87f3250() {
   return (neuron0x8d11a60()*0.695954);
}

double TMVAClassification_TMlpANN::synapse0x87f3290() {
   return (neuron0x8d11da0()*-0.282973);
}

double TMVAClassification_TMlpANN::synapse0x84a05b0() {
   return (neuron0x8066360()*-3.96987);
}

double TMVAClassification_TMlpANN::synapse0x84a05f0() {
   return (neuron0x8066790()*0.973076);
}

double TMVAClassification_TMlpANN::synapse0x8e496f0() {
   return (neuron0x8e492f0()*-0.268369);
}

double TMVAClassification_TMlpANN::synapse0x8e49730() {
   return (neuron0x8e49800()*1.2005);
}

double TMVAClassification_TMlpANN::synapse0x8e49770() {
   return (neuron0x849fee0()*0.492503);
}

double TMVAClassification_TMlpANN::synapse0x8e497b0() {
   return (neuron0x8edb310()*-0.0345909);
}

double TMVAClassification_TMlpANN::synapse0x8e40440() {
   return (neuron0x8edb7d0()*1.62609);
}

double TMVAClassification_TMlpANN::synapse0x8e40480() {
   return (neuron0x87f2dd0()*0.0599063);
}

double TMVAClassification_TMlpANN::synapse0x7b9d640() {
   return (neuron0x8066360()*2.37025);
}

double TMVAClassification_TMlpANN::synapse0x7b9d680() {
   return (neuron0x8066790()*-0.654912);
}

double TMVAClassification_TMlpANN::synapse0x84a0360() {
   return (neuron0x8e492f0()*-0.0881197);
}

double TMVAClassification_TMlpANN::synapse0x84a03a0() {
   return (neuron0x8e49800()*-0.463606);
}

double TMVAClassification_TMlpANN::synapse0x84a03e0() {
   return (neuron0x849fee0()*0.271721);
}

double TMVAClassification_TMlpANN::synapse0x84a0420() {
   return (neuron0x8edb310()*0.460572);
}

double TMVAClassification_TMlpANN::synapse0x84a0460() {
   return (neuron0x8edb7d0()*-1.35977);
}

double TMVAClassification_TMlpANN::synapse0x84a04a0() {
   return (neuron0x87f2dd0()*-0.890593);
}

double TMVAClassification_TMlpANN::synapse0x8e40bc0() {
   return (neuron0x8066360()*0.0890914);
}

double TMVAClassification_TMlpANN::synapse0x8e40c00() {
   return (neuron0x8066790()*0.200418);
}

double TMVAClassification_TMlpANN::synapse0x8e40c40() {
   return (neuron0x8e492f0()*0.133547);
}

double TMVAClassification_TMlpANN::synapse0x82dcf00() {
   return (neuron0x8e49800()*1.30143);
}

double TMVAClassification_TMlpANN::synapse0x82dcf40() {
   return (neuron0x849fee0()*0.208168);
}

double TMVAClassification_TMlpANN::synapse0x82dcf80() {
   return (neuron0x8edb310()*-0.227532);
}

double TMVAClassification_TMlpANN::synapse0x82dcfc0() {
   return (neuron0x8edb7d0()*0.859731);
}

double TMVAClassification_TMlpANN::synapse0x82dd000() {
   return (neuron0x87f2dd0()*-0.186912);
}

double TMVAClassification_TMlpANN::synapse0x8337660() {
   return (neuron0x8066360()*0.569682);
}

double TMVAClassification_TMlpANN::synapse0x83376a0() {
   return (neuron0x8066790()*-0.018496);
}

double TMVAClassification_TMlpANN::synapse0x83376e0() {
   return (neuron0x8e492f0()*-0.264944);
}

double TMVAClassification_TMlpANN::synapse0x8337720() {
   return (neuron0x8e49800()*-1.24909);
}

double TMVAClassification_TMlpANN::synapse0x8337760() {
   return (neuron0x849fee0()*0.025606);
}

double TMVAClassification_TMlpANN::synapse0x83377a0() {
   return (neuron0x8edb310()*0.19797);
}

double TMVAClassification_TMlpANN::synapse0x83377e0() {
   return (neuron0x8edb7d0()*-0.894666);
}

double TMVAClassification_TMlpANN::synapse0x8337820() {
   return (neuron0x87f2dd0()*0.0943309);
}

double TMVAClassification_TMlpANN::synapse0x8337ba0() {
   return (neuron0x8066360()*-1.87089);
}

double TMVAClassification_TMlpANN::synapse0x8337be0() {
   return (neuron0x8066790()*0.606258);
}

double TMVAClassification_TMlpANN::synapse0x8337c20() {
   return (neuron0x8e492f0()*0.255376);
}

double TMVAClassification_TMlpANN::synapse0x8337c60() {
   return (neuron0x8e49800()*1.57966);
}

double TMVAClassification_TMlpANN::synapse0x8337ca0() {
   return (neuron0x849fee0()*-1.45712);
}

double TMVAClassification_TMlpANN::synapse0x8337ce0() {
   return (neuron0x8edb310()*-0.431939);
}

double TMVAClassification_TMlpANN::synapse0x8337d20() {
   return (neuron0x8edb7d0()*1.61563);
}

double TMVAClassification_TMlpANN::synapse0x8337d60() {
   return (neuron0x87f2dd0()*1.29454);
}

double TMVAClassification_TMlpANN::synapse0x715e540() {
   return (neuron0x8066360()*-0.486779);
}

double TMVAClassification_TMlpANN::synapse0x715e580() {
   return (neuron0x8066790()*0.331817);
}

double TMVAClassification_TMlpANN::synapse0x715e5c0() {
   return (neuron0x8e492f0()*0.252336);
}

double TMVAClassification_TMlpANN::synapse0x715e600() {
   return (neuron0x8e49800()*-0.23545);
}

double TMVAClassification_TMlpANN::synapse0x715e640() {
   return (neuron0x849fee0()*-0.199996);
}

double TMVAClassification_TMlpANN::synapse0x715e680() {
   return (neuron0x8edb310()*-0.173979);
}

double TMVAClassification_TMlpANN::synapse0x715e6c0() {
   return (neuron0x8edb7d0()*0.189051);
}

double TMVAClassification_TMlpANN::synapse0x715e700() {
   return (neuron0x87f2dd0()*-0.0946229);
}

double TMVAClassification_TMlpANN::synapse0x82dd2f0() {
   return (neuron0x8066360()*-2.20847);
}

double TMVAClassification_TMlpANN::synapse0x8337da0() {
   return (neuron0x8066790()*-0.0460611);
}

double TMVAClassification_TMlpANN::synapse0x8337de0() {
   return (neuron0x8e492f0()*0.747686);
}

double TMVAClassification_TMlpANN::synapse0x87e4d20() {
   return (neuron0x8e49800()*1.07148);
}

double TMVAClassification_TMlpANN::synapse0x87e4d60() {
   return (neuron0x849fee0()*1.57593);
}

double TMVAClassification_TMlpANN::synapse0x87e4da0() {
   return (neuron0x8edb310()*-0.301566);
}

double TMVAClassification_TMlpANN::synapse0x87e4de0() {
   return (neuron0x8edb7d0()*1.32014);
}

double TMVAClassification_TMlpANN::synapse0x87e4e20() {
   return (neuron0x87f2dd0()*0.366703);
}

double TMVAClassification_TMlpANN::synapse0x87e51a0() {
   return (neuron0x87f32d0()*2.12067);
}

double TMVAClassification_TMlpANN::synapse0x87e51e0() {
   return (neuron0x8e404c0()*-0.657988);
}

double TMVAClassification_TMlpANN::synapse0x87e5220() {
   return (neuron0x8e40910()*-0.892081);
}

double TMVAClassification_TMlpANN::synapse0x87e5260() {
   return (neuron0x82dd040()*1.09359);
}

double TMVAClassification_TMlpANN::synapse0x87e52a0() {
   return (neuron0x8337860()*-1.42128);
}

double TMVAClassification_TMlpANN::synapse0x87e52e0() {
   return (neuron0x715e200()*-0.190068);
}

double TMVAClassification_TMlpANN::synapse0x87e5320() {
   return (neuron0x715e740()*1.93696);
}

