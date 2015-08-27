template <int vecLength=6>
void function_sort(float llr[vecLength], int pos[vecLength]){
	int i, j;

#pragma unroll
	for (i = 1; i < vecLength; i++) // template de classe
	{
		const float value    = llr[i];
		const int   position = pos[i];
		for (j = i; (j >= 1) && (value > llr[j - 1]); j--)
		{
			llr[j] = llr[j - 1];
			pos[j] = pos[j - 1];
		}
		llr[j] = value;
		pos[j] = position;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <int vecLength=6>
void function_sort(float illr[vecLength], float llr[vecLength], int ipos[vecLength], int pos[vecLength]){
	int i, j;
	for (i = 0; i < vecLength; i++)
	{
		llr[i] = illr[i];
		pos[i] = ipos[i];
	}

#pragma unroll
	for (i = 1; i < vecLength; i++) // template de classe
	{
		const float value    = llr[i];
		const int   position = pos[i];
		for (j = i; (j >= 1) && (value < llr[j - 1]); j--)
		{
			llr[j] = llr[j - 1];
			pos[j] = pos[j - 1];
		}
		llr[j] = value;
		pos[j] = position;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <int vecLength=6>
void function_sort(double illr[vecLength], double llr[vecLength], int ipos[vecLength], int pos[vecLength]){
	int i, j;
	for (i = 0; i < vecLength; i++)
	{
		llr[i] = illr[i];
		pos[i] = ipos[i];
	}

#pragma unroll
	for (i = 1; i < vecLength; i++) // template de classe
	{
		const float value    = llr[i];
		const int   position = pos[i];
		for (j = i; (j >= 1) && (value < llr[j - 1]); j--)
		{
			llr[j] = llr[j - 1];
			pos[j] = pos[j - 1];
		}
		llr[j] = value;
		pos[j] = position;
	}
	return;
}
