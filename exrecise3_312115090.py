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
        #group the books by title and for each group create a list of isbns.
        grouped_books = self.__books_data.groupby("Book-Title")
        





if __name__ == '__main__':
    md = myData('books.csv','ratings.csv','users.csv')
    md.mean_rating('To Kill a Mockingbird')

    