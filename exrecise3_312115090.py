# Alon Luboshitz 312115090
import sys
import pandas as pd
class myData:
    def __init__(self,books_path,ratings_path,users_path) -> None:
        self.__books_data = pd.read_csv(books_path,sep=";",encoding='latin-1',on_bad_lines='skip')
        self.__ratings_data = pd.read_csv(ratings_path,sep=";",encoding='latin-1',on_bad_lines='skip')
        self.__users_data = pd.read_csv(users_path,sep=";",encoding='latin-1',on_bad_lines='skip')
        self.remove_bad_lines()    
    def remove_bad_lines(self):
        self.__books_data.loc[:,"Year-Of-Publication"] = pd.to_numeric(self.__books_data.loc[:,"Year-Of-Publication"],errors='coerce')
        self.__books_data = self.__books_data.dropna()
        self.__books_data.loc[:,"Year-Of-Publication"] = self.__books_data.loc[:,"Year-Of-Publication"].astype('int32')
    def strip_colums(self):
       self.__ratings_data.loc[:,"ISBN"] = self.__ratings_data.loc[:,"ISBN"].str.strip().astype('str')
    def num_year(self,x,y):
        assert x < y or x ==y,('first year argument  - {} is bigger or equals to second - {}'.format(x,y))
        #assert x < 0 or y < 0,'negative year, try again'
        #filter true values of year range
        is_within_range = (self.__books_data["Year-Of-Publication"] >= x) & (self.__books_data["Year-Of-Publication"] < y)
        book_in_range = self.__books_data[is_within_range]
        print(len(book_in_range['Year-Of-Publication']))
    def df_published(self,year):
        #return spesific data frame with only books written in that year
        assert year > 0, ('negative or non year value - {}.'.format(year))
        assert isinstance(year,int),'non integer value - {}'.format(year)
        #filter only year values
        in_year = (self.__books_data['Year-Of-Publication'] == year)
        author_title = pd.DataFrame(self.__books_data[in_year], columns=['Book-Title', 'Book-Author'])
        print(author_title)
    def num_books_by_year(self,x,y):
        #filter by years
        in_range_year_df = self.__books_data[(self.__books_data["Year-Of-Publication"] >= x) & (self.__books_data["Year-Of-Publication"] < y)]
        count = in_range_year_df["Year-Of-Publication"].value_counts().sort_index()
        tuples_list = list(count.items())
        print (tuples_list)
    def mean_std(self,country):
        self.__users_data['Country'] = self.__users_data.loc[:,'Location'].str.split(',').str[-1].str.strip().astype('str')
        country_data = self.__users_data.groupby('Country').agg({'Age':['mean','std']})
        print(country_data)
        print(country_data.index)
        country_data = country_data.round(3)
        try:
            print(country_data.loc[country,'Age'])
            mean_std_tuple = tuple(country_data.loc[country,'Age'])
            return mean_std_tuple
        except KeyError:
            print("No such country: {}, exsits in data".format(country))
    def mean_rating(self,book_name):
        try:
            #group by title and get the group of the desired book
            books_title = self.__books_data.groupby('Book-Title').get_group(book_name)
        except KeyError:
            return("No book: {} have been found".format(book_name))
        #list the isbns of that book
        book_isbn = books_title["ISBN"].tolist()
        # filter the ratings data via that list
        filter_rating = self.__ratings_data[self.__ratings_data.loc[:,'ISBN'].isin(book_isbn)]
        #calculate the mean ratings after rating data been filtered.
        mean_reating = filter_rating.loc[:,'Book-Rating'].mean()
        print (mean_reating)
    '''function get k number of books rating to be presented in ascending ordrer
    thourh thier rating value'''
    def top_k(self,k):
        #check k bigger then df, k negatave etc..
        #group the books by isbn and for each group calculate mean rating.
        self.strip_colums()
        grouped_isbns = self.__ratings_data.groupby("ISBN")['Book-Rating'].mean() 
        rated_books = self.__books_data.join(grouped_isbns,on="ISBN")
        rated_books = rated_books.dropna(subset=['Book-Rating'])
        rated_books =  rated_books.sort_values(['Book-Rating','Book-Author'],ascending=[False,True])
        return rated_books.loc[:,['Book-Title','Book-Author','Book-Rating']].head(int(k))
        # pd.set_option('display.max_columns', None)
        # pd.reset_option('display.max_columns')
    '''this function group the users, sorts them by how many books they read'''
    def most_active(self,user):
        self.__ratings_data.loc[:,"User-ID"].astype('str').str.strip()
        group_id= self.__ratings_data.groupby('User-ID').size()
        group_id = group_id.sort_values(ascending=False)
        print(group_id.iloc[int(user)-1])




if __name__ == '__main__':
    md = myData('books.csv','ratings.csv','users.csv')
    md.most_active('3')

    