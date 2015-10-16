inline unsigned int LEndianReverse (unsigned int N)
{
  return ((N & 0x000000FF) << 24) |
         ((N & 0x0000FF00) << 8)  |
         ((N & 0x00FF0000) >> 8)  |
         ((N & 0xFF000000) >> 24);
}

inline unsigned short int SEndianReverse (unsigned short int N)
{
  return ((N & 0x00FF) << 8) |
         ((N & 0xFF00) >> 8);
}
