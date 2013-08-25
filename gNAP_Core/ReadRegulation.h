void ReadRegulation(double matrix[N][N])
{
	ifstream inf;
	inf.open("newGRN");
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			//fscanf(fp,"%f	",&matrix[i][j]);
			inf>>matrix[i][j];
			cout<<matrix[i][j]<<"	";
		}
		cout<<endl;
	}
	inf.close();
}