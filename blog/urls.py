from django.urls import path
from .views import (
    BlogHomeView,
    PostDetailView,
    CategoriesView,
    CategoryView,
    CreatePost,
    CreateCategoryView,
    PostDeleteView,
    PostUpdateView,
    comment_delete,
    CategoryUpdateView,
    comment_like,
    comment_dislike,
    report_comment,
    CategoryDeleteView,
    CommentDetailView,
    comment_read,
    like_post,
    dislike_post,
)

app_name = "blog"

urlpatterns = [
    path('', BlogHomeView.as_view(), name="anasayfa"),
    path('categories/', CategoriesView.as_view(), name="categories"),
    path('category/<slug:slug>/', CategoryView.as_view(), name="category"),
    path('<slug:category>/<slug:slug>/<int:pk>/<slug:author>-<int:author_id>/<created>/',
         PostDetailView.as_view(), name="post_detail"),
    path('create-post/', CreatePost.as_view(), name="post_create"),
    path('update-post/<pk>/<title>/', PostUpdateView.as_view(), name="post_update"),
    path('update-category/<pk>/<title>/', CategoryUpdateView.as_view(), name="category_update"),
    path('post_delete/<pk>/', PostDeleteView.as_view(), name="post_delete"),
    path('category_delete/<pk>/', CategoryDeleteView.as_view(), name="category_delete"),
    path('create-category/', CreateCategoryView.as_view(), name="category_create"),
    path('comment/<pk>/', CommentDetailView.as_view(), name="comment_detail"),
    path('delete/<pk>/', comment_delete, name="delete_comment"),
    path('like/<id>/', comment_like, name="comment_like"),
    path('dislike/<id>/', comment_dislike, name="comment_dislike"),
    path('report/<id>/', report_comment, name="report_comment"),
    path('okundu/<pk>/', comment_read, name="comment_read"),
    path('begendim/<pk>/<slug>/', like_post, name="like_post"),
    path('begenmedim/<pk>/<title>/', dislike_post, name="dislike_post"),

]
