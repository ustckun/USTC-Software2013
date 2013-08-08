void ReadRegulation(double matrix[170][170])
{
	ifstream inf;
	inf.open("newGRN");
	for(int i=0;i<170;i++)
	{
		for(int j=0;j<170;j++)
		{
			//fscanf(fp,"%f	",&matrix[i][j]);
			inf>>matrix[i][j];
			cout<<matrix[i][j]<<"	";
		}
		cout<<endl;
	}
	inf.close();
}