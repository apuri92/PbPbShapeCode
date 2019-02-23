double get_madness_factor(int i_cent, int i_trk, int i_jet)
{
	if (i_cent == 0)
	{
		if (i_trk == 3)
		{
			if (i_jet == 8) return 1.0027;
			if (i_jet == 9) return 1.0031;
			if (i_jet == 10) return 1.0016;
			if (i_jet == 11) return 1.0038;
			if (i_jet == 12) return 1.0020;
		}

		if (i_trk == 4)
		{
			if (i_jet == 8) return 1.0050;
			if (i_jet == 9) return 1.0051;
			if (i_jet == 10) return 1.0036;
			if (i_jet == 11) return 1.0036;
			if (i_jet == 12) return 1.0027;
		}

		if (i_trk == 5)
		{
			if (i_jet == 8) return 1.0061;
			if (i_jet == 9) return 1.0063;
			if (i_jet == 10) return 1.0023;
			if (i_jet == 11) return 1.0044;
			if (i_jet == 12) return 0.9943;
		}

		if (i_trk == 6)
		{
			if (i_jet == 8) return 1.0020;
			if (i_jet == 9) return 1.0004;
			if (i_jet == 10) return 0.9867;
			if (i_jet == 11) return 0.9777;
			if (i_jet == 12) return 0.9624;
		}

		if (i_trk == 7)
		{
			if (i_jet == 8) return 0.9968;
			if (i_jet == 9) return 0.9967;
			if (i_jet == 10) return 0.9962;
			if (i_jet == 11) return 0.9967;
			if (i_jet == 12) return 0.9992;
		}
	}
	////////
	if (i_cent == 1)
	{
		if (i_trk == 3)
		{
			if (i_jet == 8) return 1.0093;
			if (i_jet == 9) return 1.0073;
			if (i_jet == 10) return 1.0086;
			if (i_jet == 11) return 1.0088;
			if (i_jet == 12) return 1.0080;
		}

		if (i_trk == 4)
		{
			if (i_jet == 8) return 1.0114;
			if (i_jet == 9) return 1.0091;
			if (i_jet == 10) return 1.0103;
			if (i_jet == 11) return 1.0084;
			if (i_jet == 12) return 1.0144;
		}

		if (i_trk == 5)
		{
			if (i_jet == 8) return 1.0125;
			if (i_jet == 9) return 1.0111;
			if (i_jet == 10) return 1.0099;
			if (i_jet == 11) return 1.0098;
			if (i_jet == 12) return 0.9915;
		}

		if (i_trk == 6)
		{
			if (i_jet == 8) return 1.0116;
			if (i_jet == 9) return 0.9991;
			if (i_jet == 10) return 0.9917;
			if (i_jet == 11) return 0.9809;
			if (i_jet == 12) return 0.9764;
		}

		if (i_trk == 7)
		{
			if (i_jet == 8) return 0.9935;
			if (i_jet == 9) return 0.9944;
			if (i_jet == 10) return 0.9954;
			if (i_jet == 11) return 0.9942;
			if (i_jet == 12) return 0.9923;
		}
	}
	////////

	if (i_cent == 2)
	{
		if (i_trk == 3)
		{
			if (i_jet == 8) return 1.0144;
			if (i_jet == 9) return 1.0119;
			if (i_jet == 10) return 1.0111;
			if (i_jet == 11) return 1.0134;
			if (i_jet == 12) return 1.0175;
		}

		if (i_trk == 4)
		{
			if (i_jet == 8) return 1.0171;
			if (i_jet == 9) return 1.0151;
			if (i_jet == 10) return 1.0123;
			if (i_jet == 11) return 1.0168;
			if (i_jet == 12) return 1.0283;
		}

		if (i_trk == 5)
		{
			if (i_jet == 8) return 1.0223;
			if (i_jet == 9) return 1.0166;
			if (i_jet == 10) return 1.0135;
			if (i_jet == 11) return 1.0103;
			if (i_jet == 12) return 1.0237;
		}

		if (i_trk == 6)
		{
			if (i_jet == 8) return 1.0268;
			if (i_jet == 9) return 1.0076;
			if (i_jet == 10) return 1.0024;
			if (i_jet == 11) return 0.9965;
			if (i_jet == 12) return 0.9771;
		}

		if (i_trk == 7)
		{
			if (i_jet == 8) return 0.9910;
			if (i_jet == 9) return 0.9929;
			if (i_jet == 10) return 0.9929;
			if (i_jet == 11) return 0.9982;
			if (i_jet == 12) return 0.9972;
		}
	}
	////////
	if (i_cent == 3)
	{
		if (i_trk == 3)
		{
			if (i_jet == 8) return 1.0176;
			if (i_jet == 9) return 1.0177;
			if (i_jet == 10) return 1.0145;
			if (i_jet == 11) return 1.0069;
			if (i_jet == 12) return 1.0085;
		}

		if (i_trk == 4)
		{
			if (i_jet == 8) return 1.0229;
			if (i_jet == 9) return 1.0226;
			if (i_jet == 10) return 1.0192;
			if (i_jet == 11) return 1.0135;
			if (i_jet == 12) return 1.0144;
		}

		if (i_trk == 5)
		{
			if (i_jet == 8) return 1.0284;
			if (i_jet == 9) return 1.0243;
			if (i_jet == 10) return 1.0083;
			if (i_jet == 11) return 1.0079;
			if (i_jet == 12) return 1.0175;
		}

		if (i_trk == 6)
		{
			if (i_jet == 8) return 1.0467;
			if (i_jet == 9) return 1.0353;
			if (i_jet == 10) return 1.0156;
			if (i_jet == 11) return 0.9828;
			if (i_jet == 12) return 0.9884;
		}

		if (i_trk == 7)
		{
			if (i_jet == 8) return 0.9907;
			if (i_jet == 9) return 0.9908;
			if (i_jet == 10) return 0.9902;
			if (i_jet == 11) return 0.9932;
			if (i_jet == 12) return 0.9997;
		}
	}
	////////
	if (i_cent == 4)
	{
		if (i_trk == 3)
		{
			if (i_jet == 8) return 1.0229;
			if (i_jet == 9) return 1.0227;
			if (i_jet == 10) return 1.0294;
			if (i_jet == 11) return 1.0261;
			if (i_jet == 12) return 0.9960;
		}

		if (i_trk == 4)
		{
			if (i_jet == 8) return 1.0361;
			if (i_jet == 9) return 1.0341;
			if (i_jet == 10) return 1.0417;
			if (i_jet == 11) return 1.0377;
			if (i_jet == 12) return 1.0202;
		}

		if (i_trk == 5)
		{
			if (i_jet == 8) return 1.0440;
			if (i_jet == 9) return 1.0338;
			if (i_jet == 10) return  1.0390;
			if (i_jet == 11) return  1.0227;
			if (i_jet == 12) return 0.9818;
		}

		if (i_trk == 6)
		{
			if (i_jet == 8) return 1.0880;
			if (i_jet == 9) return 1.0625;
			if (i_jet == 10) return 1.0238;
			if (i_jet == 11) return 0.9894;
			if (i_jet == 12) return 0.9414;
		}

		if (i_trk == 7)
		{
			if (i_jet == 8) return 0.9919;
			if (i_jet == 9) return 0.9922;
			if (i_jet == 10) return 0.9931;
			if (i_jet == 11) return 0.9905;
			if (i_jet == 12) return 0.9913;
		}
	}
	////////
	if (i_cent == 5)
	{
		if (i_trk == 3)
		{
			if (i_jet == 8) return 1.0701;
			if (i_jet == 9) return 1.0671;
			if (i_jet == 10) return 1.0597;
			if (i_jet == 11) return 1.0438;
			if (i_jet == 12) return 1.0103;
		}

		if (i_trk == 4)
		{
			if (i_jet == 8) return 1.1077;
			if (i_jet == 9) return 1.0854;
			if (i_jet == 10) return 1.0725;
			if (i_jet == 11) return 1.0760;
			if (i_jet == 12) return 1.0174;
		}

		if (i_trk == 5)
		{
			if (i_jet == 8) return 1.1423;
			if (i_jet == 9) return 1.1085;
			if (i_jet == 10) return 1.0550;
			if (i_jet == 11) return 1.1063;
			if (i_jet == 12) return 1.0972;
		}

		if (i_trk == 6)
		{
			if (i_jet == 8) return 0.9940;
			if (i_jet == 9) return 0.9969;
			if (i_jet == 10) return 0.9955;
			if (i_jet == 11) return 0.9993;
			if (i_jet == 12) return 1.0172;
		}

		if (i_trk == 7)
		{
			if (i_jet == 8) return 0.9945;
			if (i_jet == 9) return 0.9976;
			if (i_jet == 10) return 0.9951;
			if (i_jet == 11) return 0.9998;
			if (i_jet == 12) return 1.0109;
		}
	}
	////////
}
